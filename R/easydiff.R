#' Differential Expression Analysis Using Linear Models
#'
#' This function performs differential expression analysis at the protein level using linear models.
#' It takes as input a log2-transformed expression matrix, a contrast matrix defining the comparisons of
#' interest, and group information for each sample. Optionally, a sample metadata table (\code{sample_info})
#' can be provided to specify additional information (e.g. SUBJECT identifiers). If missing, the function
#' automatically generates one assuming that each sample is from a unique biological individual.
#'
#' **Important:** The expression matrix (\code{exp_matrix}) must be log2-transformed before input. This ensures
#' that the model assumptions are met and that the computed model coefficients directly represent log2-scaled values,
#' so that the difference (fold change) between groups is computed on the log2 scale.
#'
#' The function fits a model of the form \code{lm(ABUNDANCE ~ 0 + GROUP)} for each protein, meaning that each
#' group’s coefficient corresponds directly to the mean log2 abundance in that group. Then, using the provided
#' contrast matrix, it computes the log2 fold change (logFC), standard error (SE), t-statistic (Tvalue), degrees of freedom (DF),
#' and p-value for each contrast. Missing values are omitted (no imputation is performed). In cases where a protein does not
#' have measurements in at least two groups, the function flags this with an appropriate issue message.
#'
#' @param exp_matrix A numeric matrix of expression values with proteins as rows and samples as columns.
#'   **Note:** The data must be log2-transformed. Row names should represent protein identifiers, and column names
#'   represent sample (RUN) names.
#' @param contrast.matrix A numeric matrix specifying the contrasts for comparisons. The column names must match the unique
#'   group names provided in \code{groups}. Each row of the contrast matrix defines one comparison.
#' @param groups A character vector indicating the group label for each column in \code{exp_matrix}. Its length must be equal
#'   to the number of columns in \code{exp_matrix}.
#' @param sample_info Optional. A \code{data.table} containing sample metadata. It must have at least the following columns:
#'   \itemize{
#'     \item \code{RUN}: Sample or run identifier (should match the column names of \code{exp_matrix}).
#'     \item \code{GROUP}: Group label for each sample.
#'     \item \code{SUBJECT}: Identifier for the biological sample. If technical replicates exist (multiple RUNs for the same SUBJECT),
#'       a mixed-effects model should be used with SUBJECT as a random effect. Otherwise, default is each sample is independent.
#'   }
#'   If not provided, a default \code{sample_info} is generated assuming each sample is an independent subject.
#' @param log_base Numeric. The base of the logarithm used for labeling fold change output. The default is \code{2},
#'   meaning that fold changes are reported on the log2 scale.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{\code{ComparisonResult}}{A \code{data.frame} containing the differential expression results for each protein and contrast.
#'         Columns include: \code{Protein}, \code{Label} (contrast label), fold change (named as \code{log2FC} if \code{log_base} is 2),
#'         standard error (\code{SE}), t-value (\code{Tvalue}), degrees of freedom (\code{DF}), raw p-value (\code{pvalue}),
#'         and Benjamini-Hochberg adjusted p-value (\code{adj.pvalue}). An \code{issue} column records any problems encountered
#'         (e.g. insufficient data).}
#'   \item{\code{ModelQC}}{A \code{data.frame} containing quality control metrics for each protein, including the number of
#'         measurements, missing percentage, and residuals from the fitted model. An \code{ImputationPercentage} column is also included
#'         (set to 0 since no imputation is performed).}
#'   \item{\code{FittedModel}}{A list of fitted model objects for each protein, which can be inspected further if needed.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates that \code{exp_matrix}, \code{contrast.matrix}, and \code{groups} are provided and consistent.
#'   \item If \code{sample_info} is not provided, it automatically generates sample metadata with each sample treated as an independent subject.
#'   \item Converts the expression matrix into long format and merges it with the sample metadata.
#'   \item For each protein, it fits a linear model \code{lm(ABUNDANCE ~ 0 + GROUP)}. In this model, the estimated coefficient
#'         for each group is the mean log2 abundance in that group.
#'   \item Using the provided contrast matrix, the function calculates for each contrast the log2 fold change, its standard error,
#'         t-statistic, degrees of freedom, and p-value.
#'   \item P-values are adjusted using the Benjamini-Hochberg method.
#' }
#'
#' All external functions are imported using \code{\link[=importFrom]{importFrom}}. Specifically, the following functions are used:
#' \itemize{
#'   \item \code{data.table::data.table}
#'   \item \code{data.table::melt}
#'   \item \code{data.table::merge}
#'   \item \code{data.table::rbindlist}
#'   \item \code{data.table::setnames}
#'   \item \code{data.table::as.data.table}
#'   \item \code{stats::lm}
#'   \item \code{stats::coef}
#'   \item \code{stats::vcov}
#'   \item \code{stats::residuals}
#'   \item \code{stats::fitted}
#'   \item \code{stats::pt}
#'   \item \code{stats::p.adjust}
#' }
#'
#' @examples
#' \dontrun{
#' #demo
#' set.seed(123)
#' exp_matrix <- matrix(rnorm(30), nrow = 2,
#'                      dimnames = list(c("Prot1", "Prot2"), paste0("S", 1:15)))
#'
#' log2_exp_matrix <- log2(exp_matrix + 1)
#' groups <- c(rep("health", 5), rep("diseae", 5), rep("Post", 5))
#' sample_info <- data.table::data.table(
#'   RUN = colnames(log2_exp_matrix),
#'   GROUP = groups,
#'   SUBJECT = paste0("Subject", 1:15)
#' )
#'
#' # contrast
#' contrast.matrix <- rbind(
#'   "Post/health" = c(health = -1, diseae = 0, Post = 1),
#'   "diseae/health" = c(health = -1, diseae = 1, Post = 0),
#'   "Post/diseae" = c(health = 0, diseae = -1, Post = 1)
#' )
#'
#' # diff analysis
#' result <- easydiff(log2_exp_matrix, contrast.matrix, groups, sample_info)
#'
#' head(result$ComparisonResult)# diff result
#' head(result$ModelQC)# model
#' result$FittedModel$Prot1 # protein 1
#' }
#'
#' @export
#' @importFrom data.table data.table
#' @importFrom data.table melt
#' @importFrom data.table rbindlist
#' @importFrom data.table setnames
#' @importFrom data.table as.data.table
#' @importFrom stats lm
#' @importFrom stats coef
#' @importFrom stats vcov
#' @importFrom stats residuals
#' @importFrom stats fitted
#' @importFrom stats pt
#' @importFrom stats p.adjust
#'
easydiff = function(exp_matrix,
                    contrast.matrix,
                    groups,
                    sample_info = NULL,#SUBJECT代表同一个生物学样本的多个技术重复 共享相同的SUBJECT  一般这个SUBJECT设置成不同的 开始会检查是否有重复的SUBJECT，或者是否有多个RUN属于同一个SUBJECT。如果存在这种情况，应该使用混合效应模型，其中SUBJECT作为随机效应
                    log_base = 2) {
  #不提供数据跑锤子?
  if (missing(exp_matrix) || missing(contrast.matrix) || missing(groups)) {
    stop("exp_matrix, contrast.matrix and groups are required --zcp")
  }
  if (ncol(exp_matrix) != length(groups)) {
    stop("Number of columns in exp_matrix must match length of groups --zcp")
  }
  if (!all(colnames(contrast.matrix) %in% unique(groups))) {
    stop("Column names of contrast.matrix must match group names --zcp")
  }

  #如果没有提供样本数据则自动生成 根据group和colnames(exp_matrix)一一对应生成样本数据 默认每个生物学重复只测一针即生物学样本无技术重复 SUBJECT列为不同列
  if (is.null(sample_info)) {
    sample_info = data.table::data.table(
      RUN = colnames(exp_matrix),
      GROUP = groups,
      SUBJECT = paste0("S", seq_along(groups)) #不指定就默认每个样本为独立个体
    )
  }

  #长格式
  exp_long = data.table::data.table(
    Protein = rep(rownames(exp_matrix), each = ncol(exp_matrix)),
    RUN = rep(colnames(exp_matrix), times = nrow(exp_matrix)),
    ABUNDANCE = as.vector(t(exp_matrix))
  )
  exp_long = merge(exp_long, sample_info, by = "RUN")


  model_qc_list = list()
  fitted_models = list()

  #对每个蛋白质进行差异分析
  result_list = lapply(unique(exp_long$Protein), function(prot) {
    #单个蛋白质数据
    sub_data = exp_long[Protein == prot]
    sub_data$GROUP = factor(sub_data$GROUP, levels = colnames(contrast.matrix))
    groups_present = table(sub_data$GROUP) > 0

    #记录ModelQC数据
    qc_data = sub_data[, .(
      RUN, Protein, ABUNDANCE, GROUP, SUBJECT,
      TotalGroupMeasurements = .N,
      NumMeasuredFeatures = sum(!is.na(ABUNDANCE)),
      MissingPercentage = 100 * sum(is.na(ABUNDANCE)) / .N,
      more50missing = (sum(is.na(ABUNDANCE)) / .N) > 0.5
    )]

    #检查缺失
    if (sum(groups_present) < 2) {
      res = data.table::data.table(
        Protein = prot,
        Label = rownames(contrast.matrix),
        logFC = NA_real_,
        SE = NA_real_,
        Tvalue = NA_real_,
        DF = NA_real_,
        pvalue = NA_real_,
        issue = ifelse(sum(groups_present) == 0, "allMissing",
                       "oneConditionMissing")
      )
      qc_data[, `:=`(residuals = NA, fitted = NA)]
      model_qc_list <<- c(model_qc_list, list(qc_data))
      fitted_models[[prot]] <<- NULL
      return(res)
    }

    #拟合线性模型
    model = try(lm(ABUNDANCE ~ 0 + GROUP, data = sub_data), silent = TRUE)
    if (inherits(model, "try-error")) {
      res = data.table::data.table(
        Protein = prot,
        Label = rownames(contrast.matrix),
        logFC = NA_real_,
        SE = NA_real_,
        Tvalue = NA_real_,
        DF = NA_real_,
        pvalue = NA_real_,
        issue = "modelError"
      )
      qc_data[, `:=`(residuals = NA, fitted = NA)]
      model_qc_list <<- c(model_qc_list, list(qc_data))
      fitted_models[[prot]] <<- NULL
      return(res)
    }

    #记录模型参数

    valid_rows = which(!is.na(sub_data$ABUNDANCE))
    residuals_full = rep(NA_real_, nrow(sub_data))
    residuals_full[valid_rows] = residuals(model)
    fitted_full = rep(NA_real_, nrow(sub_data))
    fitted_full[valid_rows] = fitted(model)
    qc_data[, `:=`(
      residuals = residuals_full,
      fitted = fitted_full
    )]
    model_qc_list <<- c(model_qc_list, list(qc_data))
    fitted_models[[prot]] <<- model


   # qc_data[, `:=`(
   #   residuals = residuals(model),
   #   fitted = fitted(model)
  #  )]
   # model_qc_list <<- c(model_qc_list, list(qc_data))
   # fitted_models[[prot]] <<- if (exists("model")) model else NULL

    #提取模型系数
    coefs = coef(model)
    vcov_mat = vcov(model)
    names(coefs) = gsub("^GROUP", "", names(coefs))
    colnames(vcov_mat) = rownames(vcov_mat) = gsub("^GROUP", "", colnames(vcov_mat))

    #计算每个对比
    contrast_results = lapply(1:nrow(contrast.matrix), function(i) {
      contrast_row = contrast.matrix[i, ]
      valid = names(contrast_row) %in% names(coefs)
      if (!all(valid)) {
        return(list(
          logFC = NA_real_,
          SE = NA_real_,
          Tvalue = NA_real_,
          DF = model$df.residual,
          pvalue = NA_real_,
          issue = "contrastError"
        ))
      }

      beta = coefs[names(contrast_row)]
      contrast_values = as.numeric(contrast_row)
      estimate = sum(contrast_values * beta)
      se = sqrt(t(contrast_values) %*% vcov_mat[names(contrast_row), names(contrast_row)] %*% contrast_values)
      tval = estimate / se
      pval = 2 * pt(abs(tval), df = model$df.residual, lower.tail = FALSE)

      list(
        logFC = estimate,
        SE = se,
        Tvalue = tval,
        DF = model$df.residual,
        pvalue = pval,
        issue = NA_character_
      )
    })

    #整理结果
    res = data.table::rbindlist(lapply(contrast_results, as.data.frame))
    res[, Protein := prot]
    res[, Label := rownames(contrast.matrix)]
    res[, adj.pvalue := p.adjust(pvalue, method = "BH")]
    res
  })

  #合并结果
  final_result = data.table::rbindlist(result_list, fill = TRUE)
  logFC_name = paste0("log", log_base, "FC")
  data.table::setnames(final_result, "logFC", logFC_name)

  #构建ModelQC
  model_qc = data.table::rbindlist(model_qc_list, fill = TRUE)
  model_qc[, ImputationPercentage := 0] # 假设无插补

  #包裹成List返回
  list(
    ComparisonResult = as.data.frame(final_result),
    ModelQC = as.data.frame(model_qc),
    FittedModel = fitted_models
  )
}
