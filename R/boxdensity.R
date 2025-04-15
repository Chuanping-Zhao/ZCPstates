#' Plot Boxplot and Density Plot for Quantitative QC
#'
#' This function creates a vertically stacked plot composed of a density plot and a boxplot
#' to visualize distribution and variability of intensity (e.g. protein abundance) across replicates.
#'
#' @param dt A data.frame or tibble containing the intensity values and sample information.
#' @param intensity.col Character. Name of the column representing intensity values. Default is `"Abundance"`.
#' @param sample.col Character. Name of the column representing sample or replicate labels. Default is `"BioReplicate"`.
#' @param box.outlier.color Character. Color of outliers in boxplot. Default is `"#70AD47"`.
#' @param box.fill Character. Fill color of the boxplot. Default is `"#00B050"`.
#' @param density.fill Character. Fill color of the density plot. Default is `"#00B050"`.
#' @param logTrans Logical. Whether to apply log2 transformation to the intensity column. If `FALSE`, log2 will be applied. Default is `TRUE`.
#'
#' @return A `ggplot` patchwork object combining a density plot and a boxplot.
#' 
#' @examples
#' # Example data
#' dt <- data.frame(
#'   Protein = paste0("P", 1:20),
#'   BioReplicate = rep("QC1", 20),
#'   Abundance = runif(20, 1, 1000)
#' )
#' boxdensity(dt, logTrans = FALSE)
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot stat_density theme_bw labs
#' @importFrom ggplot2 theme element_blank element_rect element_line element_text
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom rlang sym
#'
#' @export
#' 
boxdensity=function(dt,
                    intensity.col="Abundance",
                    sample.col="BioReplicate",
                    box.outlier.color="#70AD47",
                    box.fill="#00B050",
                    density.fill="#00B050",
                    logTrans=TRUE
                    ){
  
 if(!logTrans){
  dt[[intensity.col]]=log2(dt[[intensity.col]])
   }
  
  mytheme=ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),  #去除主网格线
    panel.grid.minor = ggplot2::element_blank(),  #去除次网格线
    panel.border = ggplot2::element_rect(color = "black"),  #保留边框
    axis.line = ggplot2::element_line(color = "black"),      #添加轴线
    axis.title.x = ggplot2::element_text(size = 8),  #x轴标题字体大小
    axis.title.y = ggplot2::element_text(size = 10),  #y轴标题字体大小
    axis.text.x = ggplot2::element_text(size = 8),   #x轴刻度标签字体大小
    axis.text.y =ggplot2:: element_text(size = 8) 
  )
  
  plt_qcbox=ggplot2::ggplot(data =dt ,mapping = ggplot2::aes(x=!!rlang::sym(intensity.col),y=!!rlang::sym(sample.col)))+
    ggplot2::geom_boxplot(outlier.size=0.1,
                 outlier.color=box.outlier.color,
                 show.legend = FALSE,fill=box.fill)+
    ggplot2::theme_bw()+
    ggplot2::labs(x=paste0("logTrans ",intensity.col),y="QC")+
    mytheme
    
  
  plt_density=ggplot2::ggplot(data =dt ,mapping = ggplot2::aes(x=!!rlang::sym(intensity.col)))+
    ggplot2:: stat_density(fill = density.fill,geom = "area",trim = TRUE, bw = "bcv",color="black")+
    ggplot2::labs(x=intensity.col,y="Density")+
    ggplot2::theme_bw()+
    mytheme+
    ggplot2::scale_y_continuous(expand = c(0, 0))
  

  plt=patchwork::wrap_plots(plt_density, plt_qcbox, ncol = 1, heights = c(1, 3));plt
  return(plt)  
}








