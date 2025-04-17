#' Perform differential expression analysis on a wide-format protein matrix
#'
#' This function performs linear model-based group comparison on log2-intensity
#' proteomics data. It automatically reshapes the input data, merges sample annotations,
#' handles missing value filtering, and fits fixed-effect linear models for each protein.
#'
#' @param data A wide-format data.frame or tibble of protein expression values. Rows are proteins, columns are sample Runs. One column must contain protein IDs.
#' @param contrast.matrix A contrast matrix specifying comparisons to make between experimental groups. Columns must match group names in `groupInfo$Condition`.
#' @param groupInfo A data.frame or tibble with columns `Run`, `Condition`, and `BioReplicate`, specifying the sample group assignments.
#' @param NAfilter.cutoff Optional. If specified, proteins with missing percentage above this threshold will be filtered out. (e.g. 0.5)
#' @param save_fitted_models Logical. If `TRUE`, optionally retain fitted model objects per protein. Default is `TRUE`.
#' @param protein.col Column name containing protein IDs in the input `data`. Default is `"x"`.
#' @param log2Trans Logical. If `FALSE`, input data will be log2-transformed. Default is `FALSE`.
#'
#' @return A data.table of differential expression results across all proteins, including:
#' - `Protein`: Protein ID
#' - `Label`: Contrast label (e.g. "A_vs_B")
#' - `log2FC`: Estimated log2 fold change from linear model
#' - `SE`: Standard error
#' - `Tvalue`: t statistic
#' - `DF`: Residual degrees of freedom
#' - `pvalue`: Raw p-value
#' - `adj.pvalue`: Multiple testing adjusted p-value (BH method)
#'
#' @details
#' The function uses `lm(ABUNDANCE ~ GROUP)` as the base model. Missing values are optionally filtered. Contrasts are specified manually.
#' Internally reshapes the data, merges group info, and splits by protein for independent model fitting.
#'@examples
#'# example
#'library(ZCPstates)
#'protein_wide <- data.frame(
#'x = c("P001", "P002", "P003"),
#'Sample1 = c(12.3, 13.1, 14.0),
#'Sample2 = c(12.5, 13.0, 14.2),
#'Sample3 = c(11.8, 12.9, 13.9),
#'Sample4 = c(14.1, 13.5, 15.0),
#'Sample5 = c(13.9, 13.7, 14.9),
#'Sample6 = c(14.2, 13.6, 15.1)
#')
#'sampleInfo <- data.frame(
#'Run = paste0("Sample", 1:6),
#'Condition = c("A", "A", "A", "B", "B", "B"),
#'BioReplicate = paste0("R", 1:6)
#')
#'comparison=generate_comparison_matrix(sampleInfo,
#'group_col = "Condition",
#'force_denominator = NULL)
#'
#' result <- easydiff2(
#'   data = protein_wide,
#'   contrast.matrix = comparison,
#'   groupInfo = sampleInfo,
#'   NAfilter.cutoff = 0.5,
#'   save_fitted_models = TRUE,
#'   protein.col = "x",
#'   log2Trans = TRUE
#' )
#'head(result)
#' @import data.table
#' @importFrom purrr map_dfr map_lgl
#' @importFrom stats lm coef vcov pt p.adjust model.matrix
#' @importFrom utils capture.output
#' @importFrom crayon yellow
#' @export
#'
easydiff2=function(
                          data=protein_wide,
                          contrast.matrix=comparison,
                          groupInfo=sampleInfo,
                          NAfilter.cutoff=NULL,
                         save_fitted_models = TRUE,
                         protein.col="x",
                         log2Trans=FALSE
){

  if (!protein.col %in% colnames(data)) {
    stop(paste0("The specified `protein.col` (\"", protein.col, "\") is not a valid column name in the input data.--zcp\n"))
  }


  #=====检查输入的groupInfo信息是否合法，检查输入的groupInfo信息与data中样本信息是否匹配====
  .checkGroupComparisonInput <- function(groupInfo, data) {
    required_cols <- c("Run", "Condition", "BioReplicate")

    # 检查是否包含必需列
    missing_cols <- setdiff(toupper(required_cols), toupper(colnames(groupInfo)))
    if (length(missing_cols) > 0) {
      msg <- paste0("The `groupInfo` input must include the following columns: ",
                    paste(required_cols, collapse = ", "),
                    ". Missing: ", paste(missing_cols, collapse = ", "), ".--zcp")
      stop(msg)
    }

    # 检查 groupInfo$Run 是否在 data 的列名中
    unmatched_runs <- setdiff(groupInfo$Run, colnames(data))
    if (length(unmatched_runs) > 0) {
      msg <- paste0("The following Run(s) in `groupInfo$Run` do not match any column in `data`: ",
                    paste(unmatched_runs, collapse = ", "), ".--zcp")
      stop(msg)
    }

  }

  .checkGroupComparisonInput(groupInfo,data)





  PrepareForGroupComparison=function(summarization_output,idcols=protein.col,logTrans=log2Trans){
    #summarization_output=list(protein=data, contrast.matrix=contrast.matrix,  groupInfo=groupInfo)
    #定义动态宽转长的函数
    .reshape_protein_to_long= function(proteindata_piv, protein.cols = protein.col) {
      #protein.col设为id.vars
      id_cols = protein.cols

      #宽转长
      long_dt = data.table::melt(
        data = proteindata_piv,
        id.vars = id_cols,
        variable.name = "originalRUN",
        value.name = "LogIntensities"
      )

      return(long_dt)
    }


    proteindata=data.table::as.data.table(summarization_output$protein)
    proteindata_piv=.reshape_protein_to_long(proteindata, protein.cols =protein.col)

    if (idcols != "Protein") {
      if (!idcols %in% colnames(proteindata_piv)) {
        stop(paste0("Column '", idcols, "' not found in the data."))
      }
      data.table::setnames(proteindata_piv, old = idcols, new = "Protein")
      idcols <- "Protein"  #更新变量
    }

    if (isFALSE(logTrans)) {
      proteindata_piv$LogIntensities = log2(proteindata_piv$LogIntensities)
    }



    #增加GROUP和SUBJECT信息 定义的技术重复是 GROUP_SUJECT 相同的Runs属于同一个生物分组的技术重复或者说叫生物学重复
    Info=data.table::as.data.table(summarization_output$groupInfo)
    data.table::setnames(Info,
                         old = c("Run", "Condition", "BioReplicate"),
                         new = c("originalRUN", "GROUP", "SUBJECT"))

    #合并信息
    summarized=merge(proteindata_piv, Info, by = "originalRUN", all.x = TRUE)
    summarized[, RUN := .GRP, by = originalRUN]
    #统计每个蛋白在每组数据中中的缺失值 或者全局缺失值 然后可以作为筛选阈值
    summarized[#保证NA的唯一性
      , LogIntensities := fifelse(
        LogIntensities == 0 | is.nan(LogIntensities),
        NA_real_,
        LogIntensities
      )
    ]
    summarized[#统计全局NA占比
      , MissingPercentage := mean(is.na(LogIntensities)),
      by = Protein
    ]
    summarized[#添加NA占比超过50%的蛋白信息
      , more50missing := MissingPercentage > 0.5
    ]

    group_missing = summarized[#按分组统计NA
      , .(
        GroupMissingPercentage = mean(is.na(LogIntensities)),
        GroupMore50Missing = mean(is.na(LogIntensities)) > 0.5
      ),
      by = .(Protein, GROUP)
    ]
    #合并
    summarized =merge(summarized, group_missing, by = c("Protein", "GROUP"), all.x = TRUE)


    output=split(summarized, summarized[, ..idcols])
    return(list(summarizatedoutput=output,ProteinLevelData=summarized))

  }

#把每个蛋白拆成一个List 方便建模
split_summarized = PrepareForGroupComparison(summarization_output=list(protein=data, contrast.matrix=contrast.matrix,  groupInfo=groupInfo),idcols=protein.col,logTrans=log2Trans)$summarizatedoutput
ProteinLevelData=PrepareForGroupComparison(summarization_output=list(protein=data, contrast.matrix=contrast.matrix,  groupInfo=groupInfo),idcols=protein.col,logTrans=log2Trans)$ProteinLevelData

#过滤NA
if(!is.null(NAfilter.cutoff)){
  keep_flags= purrr::map_lgl(split_summarized, function(x) {
    max(x$MissingPercentage, na.rm = TRUE) <= 0.5
  })

  filted_split_summarized =split_summarized[keep_flags]
  removed_proteins = names(split_summarized)[!keep_flags]

  message(
    "🧪 Protein filtering applied based on missing value cutoff (", NAfilter.cutoff, ").--zcp\n",
    "🚫 Proteins removed due to high missingness: ", length(removed_proteins),"--zcp\n",
    crayon::yellow( "✅ Proteins kept for modeling: ", sum(keep_flags), "--zcp\n")
  )


}else{
  filted_split_summarized =split_summarized
  message(crayon::yellow("ℹ️ No missing value filtering applied. All proteins included for modeling.--zcp\n"))
}



.checkRepeatedDesign = function(proteinleveldata) {#检查数据格式
  SUBJECT = GROUP = NULL

  input = as.data.table(proteinleveldata)
  subject_by_group = table(input[, list(SUBJECT, GROUP)])
  subject_appearances = apply(subject_by_group, 1, function(x) sum(x >  0))
  repeated = any(subject_appearances > 1)
  if (repeated) {#判断是时序设计还是对照设计 如果同一个SUBJECT 出现在多个组里面 比如 a1这个样本在 GROUP1 GROUP2 GROUP3都出现 那么他可能是时序的
    cat(crayon::yellow("Time course design of experiment\n"))
  } else {
    cat(crayon::yellow("Case control design of experiment\n"))
  }

  repeated
}
repeated=.checkRepeatedDesign(ProteinLevelData)




#为差异分析做准备
#single_protein=split_summarized$A0A075B6H7
.prepareSingleProteinForGC = function(single_protein) {
  ABUNDANCE = GROUP = SUBJECT = RUN = NULL

  data.table::setnames(single_protein,#列名一致
                       c("LogIntensities"),
                       c("ABUNDANCE"),
                       skip_absent = TRUE)
  single_protein = single_protein[!is.na(ABUNDANCE)]#去掉NA不参与建模
  single_protein[, GROUP := factor(GROUP)]#因子化
  single_protein[, SUBJECT := factor(SUBJECT)]
  single_protein[, RUN := factor(RUN)]
  return(single_protein)
}


.checkSingleSubject = function(input) {#RUN代表技术流程 比如TMT 一个RUN可能存在多个通道的样本 labelfree的实验中一个RUN就是一个样本
  SUBJECT = GROUP = NULL

  unique_annot = unique(input[, list(GROUP, SUBJECT)])
  subject_counts = unique_annot[, list(NumSubjects = data.table::uniqueN(SUBJECT)),
                                by = "GROUP"]
  all(subject_counts$NumSubject == 1)#返回FALSE代表每个GROUP下有多个生物学重复 可以建模
}


.checkTechReplicate = function(input) {#SUBJECT其实是生物学重复 比如labelfree中GROUP=A的组中 跑了3个生物学重复 那对应的SUBJECT就是A1 A2 A3
  GROUP = RUN = SUBJECT = NULL

  unique_annot = unique(input[, list(RUN,
                                     SUBJECT_NESTED = paste(GROUP,
                                                            SUBJECT,
                                                            sep = "."))])
  run_counts = unique_annot[, list(NumRuns = data.table::uniqueN(RUN)),
                            by = "SUBJECT_NESTED"]
  any(run_counts[["NumRuns"]] > 1)#返回FALSE代表每个GROUP_SUBJECT下只有一个RUN 也就是每个生物学重复不存在技术重复
}



#固定效应
fixmode_ttest = function(single_protein, contrast.matrix) {
  #因子化GROUP和改列名
  input =  .prepareSingleProteinForGC(single_protein)

  #固定效应模型
  full_fit =  lm(ABUNDANCE ~ GROUP, data = input,contrasts = list(GROUP = contr.treatment))
  beta = coef(full_fit)
  vcov_mat =  vcov(full_fit)
  df =  full_fit$df.residual

  #group level顺序和对比矩阵列名一致
  group_levels= levels(input$GROUP)
  contrast.matrix.order= contrast.matrix[, group_levels, drop = FALSE]

  #每个组的模型估计均值
  group_means= setNames(numeric(length(group_levels)), group_levels)
  for (grp in group_levels) {
    if (grp == group_levels[1]) {
      group_means[grp] <- beta["(Intercept)"]
    } else {
      coef_name <- paste0("GROUP", grp)
      group_means[grp] <- beta["(Intercept)"] + beta[coef_name]
    }
  }

  #批量计算所有的对比结果
  results =  lapply(rownames(contrast.matrix.order), function(contrast_name) {
    contrast_vector <- contrast.matrix.order[contrast_name, ]

    # group_means计算 logFC
    logFC =  sum(group_means[names(contrast_vector)] * contrast_vector)

    #L 向量 估算标准误
    L =  setNames(numeric(length(beta)), names(beta))
    for (grp in names(contrast_vector)) {
      if (grp == group_levels[1]) {
        L["(Intercept)"]  =   contrast_vector[grp]
      } else {
        coef_name=paste0("GROUP", grp)
        if (coef_name %in% names(L)) {
          L[coef_name] = contrast_vector[grp]
        } else {
          warning("Coefficient ", coef_name, " not found in model. Skipping.--zcp")
        }
      }
    }

    #计算标准误差、T值、P值
    SE=sqrt(t(L) %*% vcov_mat %*% L)
    Tvalue=logFC / SE
    pvalue=2 * pt(abs(Tvalue), df = df, lower.tail = FALSE)#双尾t 非靶是只需要查看偏离程度 无论上下调 如果事先就假设A组的蛋白比B组的高 这个时候选择单位

    data.table::data.table(
      Protein = unique(input$Protein),
      Label = contrast_name,
      log2FC = as.numeric(logFC),
      SE = as.numeric(SE),
      Tvalue = as.numeric(Tvalue),
      DF = df,
      pvalue = as.numeric(pvalue)
    )
  })

  final_result=data.table::rbindlist(results)
  final_result[, adj.pvalue := p.adjust(pvalue, method = "BH")]
  return(final_result)
}

all_results = purrr::map_dfr(filted_split_summarized, function(single_protein) {
  #跳过
  tryCatch(
    fixmode_ttest(single_protein, contrast.matrix),
    error = function(e) {
      warning("zcp:Error in protein: ", unique(single_protein$Protein), "\n", conditionMessage(e))
      return(NULL)  # 返回 NULL 会被 map_dfr 自动跳过
    }
  )
})


}
