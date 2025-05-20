#' @title Visualize Detected Proteins on HPP Abundance Curve
#'
#' @description
#' This function visualizes user-specified proteins (by UniProt ID) against the Human Protein Atlas (HPA) blood MS abundance dataset. 
#' It highlights detected proteins along a ranked abundance curve, with optional customization for point color, size, and loess smoothing.
#'
#' @param protein.ids A character vector of UniProt IDs to highlight, e.g., \code{c("P04196", "P43652", "P01871")}.
#' @param HPPms A `data.table` or `data.frame` containing HPP mass spectrometry data, typically from \code{getHPP(technology = "ms")$Data}.
#' @param pointcolor Color for highlighted points. Default is \code{"#CF3175"}.
#' @param pointalpha Transparency of points. Default is \code{0.5}.
#' @param pointsize Size of highlighted points. Default is \code{2}.
#' @param line.cex Line width of the loess smoothed curve. Default is \code{0.8}.
#' @param line.color Color of the smoothed trend line. Default is \code{"black"}.
#' @param line.alpha Transparency of the smoothed line. Default is \code{0.2}.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{Figure}}{A ggplot object showing the abundance distribution and detected proteins.}
#'   \item{\code{data}}{The input data with updated \code{detected} and \code{density} columns.}
#' }
#'
#' @details
#' Proteins are ranked by abundance, and user-supplied protein IDs are labeled as "Identified".
#' The function overlays vertical lines (by density) and points for identified proteins, along with a loess-fitted global trend.
#'
#' @examples
#' 
#' HPP_ms = ZCPstates::getHPP(dest_dir = getwd(), technology = "ms", verbose = TRUE)
#' plt= ZCPstates::easyHPPrank(
#'   protein.ids = c("P04196", "P43652", "P01871", "P10909", "P98161", 
#'                   "P05156", "Q06033", "P61073", "P05106", "P0DOX5", "Q8N6U8"),
#'   HPPms = HPP_ms$Data,
#'   pointcolor = "#CF3175",
#'   pointalpha = 0.5,
#'   pointsize = 2,
#'   line.cex = 0.8,
#'   line.color = "black",
#'   line.alpha = 0.2
#' )
#'plt$Figure
#' 
#'
#' @export
#' 
easyHPPrank=function(protein.ids,HPPms,
                     pointcolor="#CF3175",pointalpha=0.5,pointsize=2,
                     line.cex=0.8,line.color="black",line.alpha=0.2){
  protein.ids=unique(protein.ids)
  data.table::setDT(HPPms)
  HPPms[,detected:=data.table::fcase(
    Protein.Group%in%protein.ids,"Identified",
    default = "unidentified"
  )]
  
  HPPms[, density := {
    d <- density(order, bw = "nrd0")
    approx(x = d$x, y = d$y, xout = order, rule = 2)$y
  }, by = detected]
  
  marker_stats = HPPms[, .(count = .N), by = detected][
    , label := paste("Identified:", count)
  ]
  
  data.table::setDT(HPPms)
  data.table::setDT(marker_stats)
  
 plt= ggplot()+
    geom_point(aes(x=order,y=Abundance),data =HPPms[detected=="Identified"],size=pointsize,color=pointcolor,alpha=pointalpha,show.legend = FALSE)+
    geom_smooth(aes(x = order, y = Abundance),data =HPPms, method = "loess", span = 0.02,cex=line.cex, se = FALSE,alpha=line.alpha, color = line.color,show.legend = FALSE)+
    geom_segment(
      data = HPPms[detected=="Identified"],
      aes(x = order, xend = order, y = 0, yend = 2.5, alpha = density),
      color ="#238A8DFF",
      size = 0.2,
      show.legend = FALSE
    ) +
    scale_alpha(range = c(0.1, 1))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    
    xlab("Rank by Abandance")+
    ylab("logTrans[Conc pg/L]")+
    geom_text(data = marker_stats[detected=="Identified"], aes(x = Inf, y = Inf, label = label),
              hjust = 1.1, vjust = 1.1, size = 3.5, color = "black")+
    #labs(
    #   subtitle = "Coverage of proteins in the HPP (detected by MS) database across different technologies"#subtitle_text
    #  )+
    theme_bw()+
    theme(legend.position = "none",
          plot.subtitle = element_text(hjust = 0.5),
          panel.grid.major =element_blank(),
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 8),
          axis.title.y=element_text(size=10),
          axis.text.y = element_text(size = 8));plt
 
  return(list(Figure=plt,data=HPPms))
}
