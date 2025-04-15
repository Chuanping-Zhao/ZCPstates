#' Plot Correlation Matrix Heatmap with Annotations
#'
#' This function computes the correlation matrix between numeric columns of an input matrix or data frame,
#' and plots a customizable heatmap with correlation values labeled on tiles. It supports showing the full,
#' upper, or lower triangle of the correlation matrix and allows user-defined color schemes, font sizes,
#' and placement of legends.
#'
#' @param x A numeric matrix or data frame (e.g., protein abundances with samples in columns).
#' @param cor.method Character. Correlation method to use. Options: `"pearson"`, `"kendall"`, or `"spearman"`. Default is `"spearman"`.
#' @param Type Character. Type of matrix to display: `"Full"`, `"upper_tri"`, or `"lower_tri"`. Default is `"Full"`.
#' @param highcolor Character. Color for high correlation values (close to 1). Default is `"#492952"`.
#' @param lowcolor Character. Color for low correlation values (close to 0). Default is `"#1e4668"`.
#' @param midpoints Numeric. Midpoint value for color scale. Default is `0.5`.
#' @param lowpoints Numeric. Lower bound of color scale. Default is `0`.
#' @param highpoints Numeric. Upper bound of color scale. Default is `1`.
#' @param Title Character. Plot title. Default is `NULL`.
#' @param geom_tile_color Character. Color of the tile border. Default is `"white"`.
#' @param geom_tile_text_size Numeric. Border line size of tiles. Default is `0.3`.
#' @param fontcolor Character. Color of text labels on heatmap. Default is `"white"`.
#' @param lab.x.fontsize Numeric. Font size for x-axis labels. Default is `6`.
#' @param lab.y.fontsize Numeric. Font size for y-axis labels. Default is `6`.
#' @param cor.fontsize Numeric. Font size for correlation value labels. Default is `2`.
#'
#' @return A named list with:
#' \describe{
#'   \item{`plot`}{A `ggplot2` object showing the correlation heatmap.}
#'   \item{`correlation`}{The numeric correlation matrix used for plotting.}
#' }
#'
#' @examples
#' # Example
#' mat <- matrix(rnorm(100), ncol = 10)
#' colnames(mat) <- paste0("Sample", 1:10)
#' res <- heatmap_cor(mat, cor.method = "pearson", Type = "lower_tri")
#' res$plot  # Show the plot
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_text labs theme theme_classic element_text element_blank element_line element_rect scale_fill_gradient2 scale_x_discrete scale_y_discrete
#' @importFrom reshape2 melt
#' @importFrom grid unit
#' @importFrom crayon yellow
#' @export
#'
heatmap_cor <- function(x,
                        cor.method = "spearman",#cor.method = c("pearson", "kendall", "spearman")
                        Type="Full",#Type="upper_tri","lower_tri","Full"
                        highcolor= "#492952",
                        lowcolor= "#1e4668",
                        midpoints=0.5,
                        lowpoints=0,
                        highpoints=1,
                        Title=NULL,
                        geom_tile_color="white",
                        geom_tile_text_size=0.3,
                        fontcolor="white",
                        lab.x.fontsize=6,
                        lab.y.fontsize=6,
                        cor.fontsize=2
){


  if (!all(sapply(x, is.numeric))) {
    message(crayon::yellow("Warning: Non-numeric columns detected.--zcp"))
    stop("Input matrix must contain only numeric columns. --zcp")
  }


  if (is.data.frame(x)) {
    dt_cor <- x
    dt_cor[] <- lapply(dt_cor, function(x) as.numeric(as.character(x)))
  } else if (is.matrix(x)) {
    if (!is.numeric(x)) {
      stop("Input matrix must be numeric. --zcp")
    }
    dt_cor <- x
  } else {
    stop("Input must be a matrix or data.frame. --zcp")
  }


 # dt_cor=x

 # dt_cor[] <- lapply(x, function(x) as.numeric(as.character(x)))
  corvalue <- cor(as.matrix(dt_cor), method = cor.method,use = "pairwise.complete.obs")
  #提取矩阵的右下角数据函数
  get_lower_tri <- function(mat) {
    mat[upper.tri(mat)] <- NA
    return(mat)
  }
  #提取矩阵的右上角数据函数
  get_upper_tri <- function(mat) {
    mat[lower.tri(mat)] <- NA  # 将下三角部分的元素设置为NA
    return(mat)
  }
  #不对矩阵进行处理
  get_all <- function(mat){
    return(mat)
  }

  corInfo=paste0(cor.method, " correlation")


  switch (Type,
          "upper_tri" = {
            corvalue_plot <- reshape2::melt(get_upper_tri(corvalue), na.rm = TRUE) #|> filter(Var2==vars_filter)
            Figcor <- ggplot2::ggplot(corvalue_plot, ggplot2::aes(x=Var1, y=Var2, fill=value)) +
              ggplot2::geom_tile(color = geom_tile_color,size=geom_tile_text_size) +
              ggplot2::theme_classic() +
              ggplot2::labs(x = NULL,
                   y = NULL,
                   fill = corInfo,
                   title=Title) +
              ggplot2::theme(axis.text.x = ggplot2::element_text(size=lab.x.fontsize,angle = 45,vjust = 1,hjust = 0),
                    axis.text.y = ggplot2::element_text(size=lab.y.fontsize),
                    axis.title.y = ggplot2::element_text(size=8),
                    axis.title.x = ggplot2::element_text(size=8),
                    axis.ticks = ggplot2::element_blank(),
                    axis.line = ggplot2::element_blank(),
                    legend.text = ggplot2::element_text(size=8),
                    legend.title = ggplot2::element_text(size=10),
                    plot.title = ggplot2::element_text(hjust = 0.5,size=8),
                    legend.spacing.x = grid::unit(1.0, 'cm'),
                    legend.direction = "horizontal",#图例水平放置
                    legend.position = c(1, 0),     #图例放在右下角
                    legend.justification = c(1, 0)#图例左上角与1对齐
                    )+
              ggplot2::scale_fill_gradient2(high = highcolor, low = lowcolor,midpoint=midpoints,limits=c(lowpoints,highpoints), breaks = c(lowpoints, midpoints, highpoints),
                                   #labels = scales::percent_format(accuracy = 1)
              ) +
              ggplot2::geom_text(  ggplot2::aes(label = round(value, 3)),size = cor.fontsize,color=fontcolor)+
              ggplot2::scale_x_discrete(position = "top")
          },
          "lower_tri"={#获得右下角的数据
            corvalue_plot <- reshape2::melt(get_lower_tri(corvalue), na.rm = TRUE)#|> filter(Var2==vars_filter)
            Figcor <-  ggplot2::ggplot(corvalue_plot,  ggplot2::aes(x=Var1, y=Var2, fill=value)) +
              ggplot2::geom_tile(color = geom_tile_color,size=geom_tile_text_size) +
              ggplot2::theme_classic() +
              ggplot2::labs(x = NULL,
                   y = NULL,
                   fill = corInfo,
                   title=Title) +
              ggplot2:: theme(axis.text.x =  ggplot2::element_text(size=lab.x.fontsize,angle = 45,vjust = 1,hjust = 1),
                    axis.text.y = ggplot2:: element_text(size=lab.y.fontsize),
                    axis.title.y = ggplot2:: element_text(size=8),
                    axis.title.x =  ggplot2::element_text(size=8),
                    axis.ticks =  ggplot2::element_blank(),
                    axis.line = ggplot2:: element_blank(),
                    legend.text = ggplot2:: element_text(size=8),
                    legend.title = ggplot2:: element_text(size=10),
                    plot.title = ggplot2:: element_text(hjust = 0.5,size=8),
                    legend.spacing.x = grid::unit(1.0, 'cm'),
                    legend.position = "top")+
              ggplot2::scale_fill_gradient2(high = highcolor, low = lowcolor,midpoint=midpoints,limits=c(lowpoints,highpoints), breaks = c(lowpoints, midpoints, highpoints),
                                   #labels = scales::percent_format(accuracy = 1)
              ) +
              ggplot2::geom_text( ggplot2::aes(label = round(value, 2)),size = cor.fontsize,color=fontcolor)+
              ggplot2::scale_y_discrete(position = "right")
          },
          "Full"={#全部矩阵
            corvalue_plot <- reshape2::melt(get_all(corvalue), na.rm = TRUE)#|> filter(Var2==vars_filter)
            Figcor <-  ggplot2::ggplot(corvalue_plot,  ggplot2::aes(x=Var1, y=Var2, fill=value)) +
              ggplot2::geom_tile(color = geom_tile_color,size=geom_tile_text_size) +
              ggplot2::theme_classic() +
              ggplot2::labs(x = NULL,
                   y = NULL,
                   fill = corInfo,
                   title=Title) +
              theme(axis.text.x = element_text(size=lab.x.fontsize,angle = 45,vjust = 1,hjust = 1),
                    axis.text.y = element_text(size=lab.y.fontsize),
                    axis.title.y = element_text(size=8),
                    axis.title.x = element_text(size=8),
                    axis.ticks = element_blank(),
                    axis.line = element_blank(),
                    legend.text = element_text(size=8),
                    legend.title = element_text(size=10),
                    plot.title = element_text(hjust = 0.5,size=10),
                    legend.spacing.x = unit(1.0, 'cm'),
                    legend.position = "top")+
              scale_fill_gradient2(high = highcolor, low = lowcolor,midpoint=midpoints,limits=c(lowpoints,highpoints), breaks = c(lowpoints, midpoints, highpoints),
                                   #labels = scales::percent_format(accuracy = 1)
              ) +
              geom_text(aes(label = round(value, 2)),size = cor.fontsize,color=fontcolor)
          }
  )

  return(list(plot=Figcor,correlation=corvalue))

}
