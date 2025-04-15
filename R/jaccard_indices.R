#' Compute and Visualize Pairwise Jaccard Indices Between Sample Groups
#'
#' This function computes the pairwise Jaccard similarity indices of protein presence
#' across different groups (e.g., biological replicates), and visualizes the lower triangle
#' of the resulting Jaccard matrix as a heatmap.
#'
#' It removes entries from the input data where the intensity column is missing, zero,
#' "0", empty string, or NaN before computation.
#'
#' @param dt A data.frame or tibble containing at least protein identifiers, group labels, and intensity values.
#' @param intensity.col Character. Column name that stores intensity or abundance values. Default is `"Abundance"`.
#' @param protein.col Character. Column name that stores protein identifiers. Default is `"Protein"`.
#' @param group.col Character. Column name that stores sample group or replicate names. Default is `"BioReplicate"`.
#' @param fillcolor_high Character. Color used for high Jaccard values in heatmap. Default is `"#EC7607"`.
#' @param fillcolor_low Character. Color used for low Jaccard values in heatmap. Default is `"#077DEC"`.
#' @param fontsize Numeric. Font size for the Jaccard value labels shown on the heatmap tiles. Default is `2`.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{`plot`}{A ggplot object showing the lower triangle of the Jaccard index matrix.}
#'   \item{`jaccard_indices`}{A data.frame in long format representing the lower triangle of the Jaccard matrix (columns: Var1, Var2, value).}
#' }
#'
#' @examples
#' dt <- data.frame(
#'   Protein = paste0("P", 1:10),
#'   BioReplicate = rep(c("Rep1", "Rep2"), each = 5),
#'   Abundance = c(runif(5, 5, 10), runif(5, 0, 10))
#' )
#' res <- jaccard_indices(dt,
#' intensity.col = "Abundance",
#' protein.col="Protein",
#' group.col="BioReplicate",
#' fillcolor_high="#EC7607",
#' fillcolor_low="#077DEC",
#' fontsize=2)
#' res$plot  #View the heatmap
#'
#' @importFrom dplyr filter rename
#' @importFrom rlang sym
#' @importFrom ggplot2 ggplot aes geom_tile scale_x_discrete scale_y_discrete
#' @importFrom ggplot2 theme_classic theme labs element_text element_line element_blank
#' @importFrom ggplot2 scale_fill_gradient2 geom_text
#' @importFrom grid unit
#' @importFrom scales percent_format
#'
#' @export
#'
jaccard_indices=function(dt,
                         intensity.col = "Abundance",
                         protein.col="Protein",
                         group.col="BioReplicate",
                         fillcolor_high="#EC7607",
                         fillcolor_low="#077DEC",
                         fontsize=2
){
  colsym = rlang::sym(intensity.col)

 if(!is.null(intensity.col)){
   dt_Jaccard = dt |>
     dplyr::filter(!is.na(!!colsym)) |>
     dplyr::filter({
       val = !!colsym
       if (is.numeric(val)) {
         !is.nan(val) & val != 0
       } else {
         val != "" & val != "0"
       }
     })




 }else{
   dt_Jaccard = dt
 }


  Jaccard_set_list=split(dt_Jaccard[[protein.col]], dt_Jaccard[[group.col]])


  .calculate_jaccard =function(set1, set2) {
    intersection = length(intersect(set1, set2))
    union = length(union(set1, set2))
    jaccard_index =intersection / union
    return(jaccard_index)
  }

  jaccard_indices= matrix(nrow = length(Jaccard_set_list), ncol = length(Jaccard_set_list))
  #计算每一对的 Jaccard 指数
  for (i in seq_along(Jaccard_set_list)) {
    for (j in seq_along(Jaccard_set_list)) {
      jaccard_indices[i, j] = .calculate_jaccard(Jaccard_set_list[[i]], Jaccard_set_list[[j]])
    }
  }

  # 给矩阵加上行名和列名
  rownames(jaccard_indices)= colnames(jaccard_indices)= names(Jaccard_set_list)

  #提取矩阵的右下角
  get_lower_tri = function(mat) {
    mat[upper.tri(mat)] = NA
    return(mat)
  }
  #获得右下角矩阵数据
  #jaccard_indices_plot = reshape2::melt(get_lower_tri(jaccard_indices), na.rm = TRUE)


  jaccard_indices_plot= get_lower_tri(jaccard_indices) |>
    as.table() |>
    as.data.frame() |>
    dplyr::filter(!is.na(Freq)) |>
    dplyr::rename(Var1 = Var1, Var2 = Var2, value = Freq)





  #定义标签
  labz =unique(jaccard_indices_plot$Var1)

  jaccard_indices_mean =mean(jaccard_indices_plot$value)



  plt_jaccard_indices = ggplot2::ggplot(jaccard_indices_plot,ggplot2:: aes(x=Var1, y=Var2, fill=value)) +
    ggplot2::geom_tile(color = "black",size=0.3)+
    ggplot2::scale_x_discrete(labels= labz) +
    ggplot2::scale_y_discrete(labels= labz, position = "right") + #y轴置于右侧
    ggplot2::theme_classic() +
    ggplot2::labs(x =NULL,
         y = NULL,
         fill = "Jaccard Index",
         title="Protein Jaccard Index",
         subtitle = paste0("Average Jaccard Index:",round(jaccard_indices_mean, 2)*100,"%")) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 5,angle = 45,vjust = 1,hjust = 1),
          axis.text.y = ggplot2::element_text(size=5),
          axis.title.y = ggplot2::element_text(size=5),
          axis.title.x = ggplot2::element_text(size=5),
          legend.text = ggplot2::element_text(size=10),
          legend.title = ggplot2::element_text(size=10),
          plot.title = ggplot2::element_text(hjust = 0.5,size=10),
          plot.subtitle = ggplot2::element_text(hjust = 0.5,size=10,color = "#EC7607"),
          axis.line = ggplot2::element_line(size = 0.3),  # 设置坐标轴线的粗细
          axis.ticks = ggplot2::element_line(size = 0.3),  # 设置刻度线的粗细
          legend.spacing.x = grid::unit(1.0, 'cm'),
          legend.position = "top")+
    ggplot2:: scale_fill_gradient2(high = fillcolor_high, low =fillcolor_low ,midpoint=0.5,limits=c(0, 1), breaks = c(0, 0.5, 1),
                         labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::geom_text(ggplot2::aes(label = round(value, 2)*100),size = fontsize);plt_jaccard_indices

  return(list(
    plot = plt_jaccard_indices,
    jaccard_indices = jaccard_indices_plot
  ))


  }
