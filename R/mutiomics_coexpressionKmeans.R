#' Multi-Omics Coexpression Clustering and Visualization with KMeans
#'
#' This function performs multi-omics coexpression analysis between two omics datasets
#' (e.g., metabolome and proteome), including group-wise mean calculation, KMeans clustering,
#' pairwise correlation, and aligned cluster visualization.
#'
#' @param omics1_exp A data frame of the upstream omics expression matrix (e.g., metabolomics). Rows are features (e.g., metabolites), columns are samples.
#' @param omics2_exp A data frame of the downstream omics expression matrix (e.g., proteomics). Rows are features (e.g., proteins), columns are samples.
#' @param omics1_type A character string representing the name of the upstream omics type. Default is `"Metabolome"`.
#' @param omics2_type A character string representing the name of the downstream omics type. Default is `"Proteome"`.
#' @param groupInfo A data frame with exactly **two columns**: `Condition` and `sample`. `Condition` is the group/condition name; `sample` corresponds to column names of `omics1_exp` and `omics2_exp`.
#' @param psych.corr.test.method Correlation method to use (`"pearson"` or `"spearman"`). Default is `"pearson"`.
#' @param pvaluecutoff Significance cutoff for correlation p-value. Default is `0.01`.
#' @param coefficientcutoff Absolute correlation coefficient cutoff. Default is `0.8`.
#' @param clustercounts Number of clusters to use in KMeans. If `NULL`, the optimal number will be estimated automatically via `cascadeKM`. Default is `NULL`.
#' @param positive_only Logical. If `TRUE`, only proteins positively correlated with metabolites are retained. Default is `TRUE`.
#' @param omics1.linecolor Color for individual lines in upstream omics plots. Default is `"#ffb6ce"`.
#' @param omics1.meancolor Color for mean line in upstream omics plots. Default is `"#E1151D"`.
#' @param omics2.linecolor Color for individual lines in downstream omics plots. Default is `"#7895c3"`.
#' @param omics2.meancolor Color for mean line in downstream omics plots. Default is `"#103883"`.
#'
#' @return A list containing:
#' \item{plot}{A patchwork plot object visualizing aligned clusters of omics1 and omics2.}
#' \item{omics1cluster}{Clustered upstream omics (omics1) data with z-score normalization.}
#' \item{omics2_coexpression_cluster}{Filtered and aligned downstream omics (omics2) data based on coexpression.}
#' \item{mutiomics_coexpression}{A data frame of significant coexpression pairs and their correlation statistics.}
#'
#' @importFrom psych corr.test
#' @importFrom vegan cascadeKM
#' @importFrom dplyr left_join arrange distinct
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line theme element_blank element_text element_line element_rect
#' @importFrom cowplot get_legend plot_grid
#' @importFrom patchwork plot_layout
#' @importFrom purrr reduce
#'
#' @examples
#' \dontrun{
#' data("metaRawdata")
#' data("proteinRawdata")
#' data("sampleInfo")
#'
#' meta_expr <- metaRawdata_pos |> column_to_rownames("ID")
#' protein_expr <- proteinRawdata |> column_to_rownames("ID")
#'
#' result <- mutiomics_coexpressionKmeans(
#'   omics1_exp = meta_expr,
#'   omics2_exp = protein_expr,
#'   groupInfo = sampleInfo,
#'   coefficientcutoff = 0.9,
#'   clustercounts = NULL,
#'   positive_only = TRUE
#' )
#'
#' # View plot
#' result$plot
#' # Save result
#' saveRDS(result, file = "coexpression_result.rds")
#' }
#'
#' @export
#'
mutiomics_coexpressionKmeans=function(
    omics1_exp, # omics1_exp=metaRawdata_pos |> column_to_rownames(var = "ID")
    omics2_exp, # omics2_exp=proteinRawdata |> column_to_rownames(var = "ID")
    omics1_type="Metabolome",
    omics2_type="Proteome",
    groupInfo,  # groupInfo=sampleInfo
    psych.corr.test.method="pearson",
    pvaluecutoff=0.01,
    coefficientcutoff=0.8,
    clustercounts=NULL,#默认NULL代表自动选择最佳聚类数量 也可以人为指定
    positive_only = TRUE,#仅保存与差异代谢物变化趋势正相关的蛋白
    omics1.linecolor= '#ffb6ce',
    omics1.meancolor = '#E1151D',

    omics2.linecolor = '#7895c3',
    omics2.meancolor = '#103883'
    ){
  # omics1_exp代谢组的原始矩阵 没有做logTrans
  # omics2_exp蛋白质组的原始矩阵 没有做logTrans


  options(scipen = 9, stringsAsFactors = F)
  omics1_exp[is.na(omics1_exp)]=0
  omics2_exp[is.na(omics2_exp)]=0


  #  quantFrame=omics1_exp
  #  sampleFrame=groupInfo

  .calcMeans = function(quantFrame, sampleFrame) {
    groupNames =unique(sampleFrame$Condition)
    for (group in groupNames) {
      samples = sampleFrame[sampleFrame$Condition == group, ]$sample
      #quantFrame= lapply(quantFrame, as.numeric)

      quantFrame[, samples] = lapply(quantFrame[, samples], as.numeric)


      quantFrame$group = apply(quantFrame[, as.character(samples), drop = FALSE],
                               1, mean,na.rm = TRUE)
      names(quantFrame)[names(quantFrame) == "group"] <- group
    }
    # Select neccessary coloums
    subQuant = subset(quantFrame, drop = F, select = as.character(groupNames))
    subQuant = as.data.frame(subQuant, check.names = FALSE)
    return(subQuant)
  }

  omics1GroupDF=.calcMeans(omics1_exp, groupInfo)#--代谢组学均值
  omics2GroupDF=.calcMeans(omics2_exp, groupInfo)#--蛋白组学均值

  # metaclusterFrame <- kmeansAnalysis(metaGroupDF, 10)
  # quantFrame=omics2GroupDF;clusterNum=10
  kmeansAnalysis = function(quantFrame,clusters) {

    .scaleQuant <- function(quantFrame) {
      quantTranspose <- t(as.matrix(quantFrame))

      quantConvert <- quantTranspose  |>
        scale(center = T, scale = T) |>
        t()|>
        na.omit()
      return (quantConvert)
    }
    quantConvert = .scaleQuant(log2(quantFrame))
    # calculate the best number of clusters
    set.seed(20250424)


    if(is.null(clusters)){
      kmean.fit =vegan::cascadeKM(
        data = quantConvert,
        inf.gr = 1,  # The number of groups for the partition with the smallest number of groups of the cascade (min).
        sup.gr = 10,  # The number of groups for the partition with the largest number of groups of the cascade (max).
        iter = 100,#The number of random starting configurations for each value of K
        criterion= "calinski")  # ""ssi" or "calinski"
      plot(kmean.fit)
      clusterNum= kmean.fit$results[2, ]  |>   which.max()  |>   as.numeric()
    }else{
      clusterNum=clusters
    }




    ##　clusterNum为K值，可以自行定义，也可以使用前面老师讲解的方法计算最佳K值
    quantCluster =kmeans(quantConvert, centers = clusterNum)
    ## 将Kmeans分析后的代谢物对应的聚类簇与z-score数据合并成数据框
    clusterFrame = data.frame(`cluster` = quantCluster$cluster,
                              quantConvert, check.names = FALSE) |>
      dplyr::arrange(cluster) |>  # 排序
      tibble::rownames_to_column(var = 'Index') # 行名转化为Index列
    return(clusterFrame)
  }


  omics1clusterFrame = kmeansAnalysis(omics1GroupDF,clusters=clustercounts)#-上游组学kmean
  omics2clusterFrame = kmeansAnalysis(omics2GroupDF,clusters=clustercounts)#-下游组学keman 未筛选

  # 上游组学Zscore
  ## scale
  .scaleQuant <- function(quantFrame) {
    quantTranspose <- t(as.matrix(quantFrame))

    quantConvert <- quantTranspose  |>
      scale(center = TRUE, scale = TRUE) |>
      t()|>
      na.omit()
    return (quantConvert)
  }
  proteinFrame =data.frame(.scaleQuant(omics2GroupDF), check.names = FALSE)  |>
    tibble::rownames_to_column(var = 'Index')

  # 相关性
  omics1_exp_numeric = as.data.frame(lapply(omics1_exp, as.numeric))
  omics2_exp_numeric =as.data.frame(lapply(omics2_exp, as.numeric))
  #加上行名
  omics1_rownames=rownames(omics1_exp)
  omics2_rownames=rownames(omics2_exp)
  rownames(omics1_exp_numeric)=omics1_rownames
  rownames(omics2_exp_numeric)= omics2_rownames

  corres=psych::corr.test(t(omics1_exp_numeric), t(omics2_exp_numeric), method = psych.corr.test.method, adjust = 'none', ci = FALSE)
  corCmat=corres$r
  corPmat = corres$p
  ## 合并相关系数R和Pvalue，即使每行表示一个代谢物与某个基因的相关系数和P值
  res =data.frame(Index1 = rep(rownames(corCmat), ncol(corCmat)),
                  Index2 = rep(colnames(corCmat), each = nrow(corCmat)),
                  coefficient = c(corCmat), pvalue = c(corPmat))
  #res= res[abs(res[, 3]) >= coefficientcutoff & res[, 4] <= pvaluecutoff, ,drop = FALSE] # 筛选相关性系数大于0.8且P<0.05的蛋白-代谢物
 # res= res[res[, 3] >= coefficientcutoff & res[, 4] <= pvaluecutoff, ,drop = FALSE]  #只筛选正相关

  if (positive_only) {
    res = res[res$coefficient >= coefficientcutoff & res$pvalue <= pvaluecutoff, ,drop = FALSE]
  } else {
    res = res[abs(res$coefficient) >= coefficientcutoff & res$pvalue <= pvaluecutoff, ,drop = FALSE]
  }



  ## 获取上游组学的聚类簇
  cogene = res  |>
    dplyr::left_join(subset(omics1clusterFrame, select = c('Index', 'cluster')),
                     by = c('Index1' = 'Index'))  |> # 合并代谢组[下游组学及其对应的上游组学共表达显著的蛋白]
    dplyr::arrange(cluster, dplyr::desc(abs(coefficient)))  |>  # 排序，聚类簇升序，相关性绝对值倒序
    dplyr::distinct(Index2, cluster, .keep_all = FALSE)

  ## 获取上游组学蛋白的Z-Score
  geneclusterFrame =cogene  |>  dplyr::left_join(proteinFrame, by = c('Index2' = 'Index'))  |>  # 合并 与代谢共表达显著的蛋白信息 与 蛋白均值矩阵
    dplyr::rename(Index = Index2)






  #绘图
  plotplot =function(df, clusterNum, show_x = TRUE, show_y = TRUE,
                     linecolor = 'red', meancolor = 'orange') {
    groups = colnames(df)[3:ncol(df)]
    anadata = df[df$cluster == clusterNum, c(1, 3:ncol(df))]
    plotdata = tidyr::pivot_longer(anadata, !Index, names_to = "Group", values_to = "value")
    plotdata$Group = factor(plotdata$Group, levels = groups, ordered = TRUE)
    meandata = data.frame(Group = groups, value = colMeans(anadata[, groups]))
    meandata$Group = factor(meandata$Group, levels = groups, ordered = TRUE)

    p = ggplot2::ggplot(plotdata, ggplot2::aes(x = Group, y = value, group = Index)) +
      ggplot2::geom_line(color = linecolor) +
      ggplot2::geom_line(data = meandata, ggplot2::aes(x = Group, y = value, group = 1),
                         color = meancolor, linewidth = 1) +
      ggplot2::labs(title = paste0("Cluster ", clusterNum, " (", nrow(anadata), ")")) +
      ggplot2::xlab(NULL) +
      ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(color = 'black', linewidth = 0.5),
                     plot.title = ggplot2::element_text(face = 'plain', hjust = 0.5,size=10),
                     panel.background = ggplot2::element_rect(fill = 'white'))

    if (!show_x) {
      p = p + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                             axis.ticks.x = ggplot2::element_blank())
    } else {
      p = p + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, angle = 65, hjust = 1))
    }

    if (!show_y) {
      p= p + ggplot2::ylab(NULL) +
        ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank())
    } else {
      p = p + ggplot2::labs(y = 'Z-score') +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8),
                       axis.title.y = ggplot2::element_text(size = 8))
    }

    return(p)
  }


  common_clusters = sort(intersect(
    unique(omics1clusterFrame$cluster),
    unique(geneclusterFrame$cluster)
  ))

  plist = list()
  n_clusters = length(common_clusters)

  for (i in seq_along(common_clusters)) {
    clusterNum = common_clusters[i]


    p1 = plotplot(omics1clusterFrame, clusterNum,
                  show_x = FALSE,
                  show_y = (i == 1),
                  linecolor = omics1.linecolor,
                  meancolor = omics1.meancolor)

    p2 = plotplot(geneclusterFrame, clusterNum,
                  show_x = TRUE,
                  show_y = (i == 1),
                  linecolor = omics2.linecolor,
                  meancolor = omics2.meancolor)

    # 每列拼为一个组
    plist[[i]] = p1 / p2
  }


  #增加图例
  legend_df =data.frame(
    Group = c(omics1_type, omics2_type),
    value = 1,
    type =  c(omics1_type, omics2_type)
  )
  legend_colors = setNames(
    c(omics1.meancolor, omics2.meancolor),
    c(omics1_type, omics2_type)
  )

  legend_plot = ggplot2::ggplot(legend_df, ggplot2::aes(x = Group, y = value, color = type)) +
    ggplot2:: geom_point(size = 5) +
    ggplot2::scale_color_manual(values = legend_colors, name = "Type") +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right",
                   legend.title =ggplot2:: element_blank(),
                   legend.text = ggplot2::element_text(size = 8))
  legend_only = cowplot::get_legend(legend_plot)


  patch_plot = purrr::reduce(plist, `|`) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "right")

  final_plot=cowplot::plot_grid(patch_plot, legend_only, ncol = 2, rel_widths = c(1, 0.15))
  final_plot

  return(
    list(
      plot=final_plot,
      omics1cluster=omics1clusterFrame,
      omics2_coexpression_cluster=geneclusterFrame,
      mutiomics_coexpression=res
    )
  )



}
