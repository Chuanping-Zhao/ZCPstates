#' Visualize easyDiff Results with a Volcano Plot
#'
#' `easyvolcano` generates volcano plots to visualize the output of differential analysis results from `easydiff`.
#' It highlights features (e.g., proteins, peptides, metabolites) that show significant changes between conditions,
#' based on specified fold-change and significance thresholds.
#'
#' @param diff The result of an `easydiff` analysis, provided as a tibble or data frame with fixed columns:
#'   \code{Label}, \code{Protein}, \code{log2FC}, and \code{pvalue}. Each row represents a feature (e.g., a protein or metabolite),
#'   with its identifier (Label), associated protein or group (Protein), log2 fold-change (log2FC), and p-value (pvalue) from the differential analysis.
#' @param cutoff_fc A numeric threshold for the fold-change cutoff (on the log2 scale). Features with an absolute log2 fold-change
#'   greater than or equal to this value are considered to have significant fold changes. Default is \code{1.2}.
#' @param cutoff_p A numeric threshold for the significance level (p-value). Features with p-values less than or equal to this
#'   threshold are considered statistically significant. Default is \code{0.05}.
#' @param top_marker An integer specifying the number of top significant features to label on the plot. The function will label
#'   the \code{top_marker} features with the smallest p-values (most significant) in the volcano plot (e.g., use 10 to label the top 10 most significant features).
#' @param max_overlaps A value controlling the maximum number of overlapping labels allowed, passed to \code{ggrepel::geom_label_repel()}.
#'   This helps avoid clutter by limiting how many labels overlap on the plot. Adjust this value to fine-tune label placement.
#' @param Feature A string indicating the type of feature being analyzed (e.g., \code{"protein"}, \code{"glycopeptide"}, \code{"metabolite"}).
#'   This is used in plot annotations (such as axis labels or titles) to clarify what type of features are being visualized.
#'
#' @importFrom dplyr filter arrange mutate
#' @importFrom ggplot2 ggplot aes geom_point labs theme scale_color_manual scale_x_continuous scale_y_continuous
#' @importFrom ggrepel geom_text_repel
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{pltvol1}}{A \code{ggplot2} object showing a standard volcano plot (log2 fold change on the x-axis
#'     and -log10 p-value on the y-axis) of the `easydiff` results, with points highlighted based on significance and the top features labeled.}
#'   \item{\code{pltvol2}}{A \code{ggplot2} object providing an alternative volcano plot visualization of the results (for example, a variation
#'     in styling or labeling compared to \code{pltvol1}). This plot also highlights significant features and labels the top features.}
#'   \item{\code{dt_vol}}{A filtered data frame (tibble) containing only the significant features—those that meet the
#'     \code{cutoff_fc} and \code{cutoff_p} thresholds. This can be used for further analysis or inspection of the significant features.}
#' }
#'
#' @examples
#' \dontrun{
#' # Suppose `diff_results` is a tibble from an easyDiff analysis:
#' diff_results <- easydiff(diff.matrix.log2trans, contrast.matrix,groups,sample_info)
#'
#' # Generate volcano plots with default thresholds and label top 10 proteins:
#' volcano_output <- easyvolcano(diff_results, cutoff_fc = 1.2, cutoff_p = 0.05, top_marker = 10, Feature = "protein")
#'
#' # The result is a list. We can display the two volcano plots:
#' print(volcano_output$pltvol1)  # Primary volcano plot
#' print(volcano_output$pltvol2)  # Alternative volcano plot
#'
#' # View the data frame of significant features:
#' head(volcano_output$dt_vol)
#' }
#'
#' @export
#'
easyvolcano=function(diff,cutoff_fc=1.2,cutoff_p=0.05,top_marker=0,max_overlaps=0,Feature="protein"){

  diff.result=diff |> dplyr::mutate(significance= dplyr::case_when(
    (log2FC>=log2(cutoff_fc))&(pvalue<=cutoff_p)~"up",
    (log2FC<=-log2(cutoff_fc))&(pvalue<=cutoff_p)~"down",
    TRUE~"nosig"
  ))


down_up_counts=diff.result  |>dplyr:: filter(significance!="nosig")|> dplyr::group_by(Label,significance) |> dplyr::summarise(cout=dplyr::n())

plt_vol=
  ggplot2::ggplot()+
  ggplot2:: geom_point(data=diff.result |> dplyr::filter(significance=="nosig"),ggplot2::aes(x=log2FC, y=-log10(pvalue)),size = 3,alpha = 0.1,shape = 21,color = "black",fill="black",stroke = 0.5)+

  ggplot2:: geom_point(data=diff.result |> dplyr::filter(significance=="up"),ggplot2::aes(x=log2FC, y=-log10(pvalue)),size = 3,alpha = 0.8,shape = 21,color = "black",fill="#A40000",stroke = 0.5)+

  ggplot2::geom_point(data=diff.result |> dplyr::filter(significance=="down"),ggplot2::aes(x=log2FC, y=-log10(pvalue)),size = 3,alpha = 0.8,shape = 21,color = "black",fill="#4E708A",stroke = 0.5)+

  ggplot2:: labs(x = "log2(Fold change)",
       y =paste0("-log10(pvalue)"),
       title =paste0("Volcano of ",Feature) ) +
  ggplot2::geom_hline(yintercept = c(-log10(cutoff_p)),
             linewidth = 0.5,
             color = "#3288BD",
             lty = "dashed")+
  ggplot2::geom_vline(xintercept = c(log2(cutoff_fc),log2(1/cutoff_fc)),
             linewidth = 0.5,
             color = "#3288BD",
             lty = "dashed")+
  ggplot2:: geom_text(data = down_up_counts |> dplyr:: filter(significance == "up"),
                      ggplot2::aes(x = Inf, y = Inf, label = paste0("n(up)=",cout)),
            vjust = 2, hjust = 1, size = 3, color = "#A40000", inherit.aes = FALSE) +
  ggplot2::geom_text(data = down_up_counts |> dplyr:: filter(significance == "down"),
                     ggplot2::aes(x = -Inf, y = Inf, label  = paste0("n(down)=",cout)),
            vjust = 2, hjust = 0, size = 3, color = "#4874CB", inherit.aes = FALSE)+

  ggplot2::facet_wrap(~Label,scales = "free")+
  ggplot2::theme_bw(base_size = 10) +
  #ggplot2::xlim(-3, 3)+
 # ggplot2::ylim(0, 10)+
  ggplot2::theme(plot.title = ggplot2::element_text(size = 15, hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5),
        plot.background = ggplot2::element_rect(fill = "white", color = NA),
        panel.background = ggplot2::element_rect(fill = "white", color = NA),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        panel.ontop = FALSE)

#vol2
dat.plot=diff.result |> dplyr::mutate(Label=factor(diff.result$Label,levels = unique(diff.result$Label)))






background.dat <- data.frame(
  dat.plot |>  dplyr::group_by(Label) |> dplyr:: filter(log2FC>0) |>
    dplyr::summarise("y.localup"=max(log2FC)),
  dat.plot  |> dplyr:: group_by(Label) |> dplyr::filter(log2FC<0) |>
    dplyr::summarise("y.localdown"=min(log2FC)),
  x.local=seq(1:length(unique(dat.plot$Label)))
) |> dplyr:: select(-Label.1)


x.number <- background.dat  |>   dplyr:: select(Label,x.local)
dat.plot <- dat.plot |>   dplyr::left_join(x.number,by = "Label")

#selecting top-up and top-down proteins
dat.marked.up <- dat.plot  |>  dplyr::filter(significance=="up")  |>
  dplyr::group_by(Label)  |> dplyr::arrange(-log2FC)  |>
  dplyr::top_n(top_marker,abs(log2FC))
dat.marked.down <- dat.plot  |>  dplyr::filter(significance=="down") |>
  dplyr::group_by(Label)  |>  dplyr::arrange(log2FC)  |>
  dplyr:: top_n(top_marker,abs(log2FC))
dat.marked <- dat.marked.up  |>  dplyr::bind_rows(dat.marked.down)
#referring group information data
dat.infor <- background.dat  |>
  dplyr::mutate("y.infor"=rep(0,length(Label)))

dat.count= dat.infor|>
  dplyr::left_join(down_up_counts,by = "Label")

color.pals=c(
  "#6181BD4E",
  "#F348004E",
  "#64A10E4E",
  "#9300264E",
  "#464E044E",
  "#049a0b4E",
  "#4E0C664E",
  "#D000004E",
  "#FF6C004E",
  "#FF00FF4E",
  "#c7475b4E",
  "#00F5FF4E",
  "#BDA5004E",
  "#A5CFED4E",
  "#f0301c4E",
  "#2B8BC34E",
  "#FDA1004E",
  "#54adf54E",
  "#CDD7E24E",
  "#9295C14E"
)




plt.vol2=
  ggplot2::ggplot()+
  ggplot2::geom_col(background.dat,mapping= ggplot2::aes(x.local,y.localup),
                    fill="grey50",alpha=0.2,width=0.9,just = 0.5)+
  ggplot2::geom_col(background.dat,mapping= ggplot2::aes(x.local,y.localdown),
                    fill="grey50",alpha=0.2,width=0.9,just = 0.5)+
  ggplot2::geom_jitter(dat.plot |> dplyr::filter(significance=="nosig"),mapping= ggplot2::aes(x.local,log2FC),color="#eaebea",size=1.5,width = 0.4,alpha= 0.8)+#color="#eaebea"
  ggplot2::geom_jitter(dat.plot |> dplyr::filter(significance=="up"),mapping= ggplot2::aes(x.local,log2FC+0.5*(1/cutoff_fc)),color="#d56e5e",size=1.5,width = 0.4,alpha= 0.4)+
  ggplot2::geom_jitter(dat.plot |> dplyr::filter(significance=="down"),mapping= ggplot2::aes(x.local,log2FC-0.5*(1/cutoff_fc)),color="#5390b5",size=1.5,width = 0.4,alpha= 0.4)+

  ggplot2::geom_tile(dat.infor,mapping= ggplot2::aes(x.local,y.infor,fill=Label,color = Label),show.legend = F,
                     height=1.5,
                     color = color.pals[1:length(unique(dat.plot$Label))],
                     fill = color.pals[1:length(unique(dat.plot$Label))],
                     alpha = 0.6,
                     width=0.9)+
  ggplot2::guides(size= ggplot2::guide_legend(title="Count"))+
  ggplot2::labs(x=NULL,y="log2 Fold change",
                subtitle = paste0("|log2FC|>=",round(log2(cutoff_fc),3)," & p value<",cutoff_p))+
  ggplot2::geom_text(dat.infor,mapping= ggplot2::aes(x.local,y.infor,label=Label),show.legend = F)+

  ggplot2:: geom_text(data = dat.count |> dplyr:: filter(significance == "up"),
                      ggplot2::aes(x = x.local, y = Inf, label = paste0("n(up)=",cout)),
                      vjust = 1, hjust = 0.5, size = 3, color = "#A40000", inherit.aes = FALSE) +
  ggplot2::geom_text(data = dat.count |> dplyr:: filter(significance == "down"),
                     ggplot2::aes(x = x.local, y = -Inf, label  = paste0("n(down)=",cout)),
                     vjust = -1, hjust =0.5  , size = 3, color = "#4874CB", inherit.aes = FALSE) +
  ggrepel::geom_text_repel(dat.marked.up,mapping= ggplot2::aes(x.local,log2FC,label=Protein),
                           color="#d56e5e",
                           #force = 0.5,
                           size=2,
                           show.legend = F,
                           max.overlaps = max_overlaps,
                           seed = 233,
                           min.segment.length = 0.2,
                           force_pull = 2,
                           box.padding = 0.1,
                           segment.linetype = 1,
                           segment.color = 'black',
                           segment.alpha = 0.5,
                           # direction = "y",
                           nudge_y = log2(cutoff_fc),
                           #hjust = 0.5
  )+
  ggrepel::geom_text_repel(dat.marked.down,mapping= ggplot2::aes(x.local,log2FC,label=Protein),
                           color="#5390b5",
                           show.legend = F,
                           # force = 0.5,
                           size=2,
                           max.overlaps = max_overlaps,
                           seed = 233,
                           min.segment.length = 0.2,
                           #force_pull = 2,
                           box.padding = 0.1,
                           segment.linetype = 1,
                           segment.color = 'black',
                           segment.alpha = 0.5,
                           nudge_y = -log2(cutoff_fc),
                           #direction = "y",
                           #hjust = 0.5
  )+
  ggplot2::theme_classic()+
  ggplot2::theme(
    axis.ticks.x =  ggplot2::element_blank(),
    axis.ticks.y = ggplot2:: element_blank(),
    axis.line.x = ggplot2:: element_blank(),
    axis.text.x = ggplot2:: element_blank(),
    axis.text.y = ggplot2:: element_blank(),
    legend.spacing.x= ggplot2::unit(0.1,'cm'),
    legend.key.width= ggplot2::unit(0.5,'cm'),
    legend.key.height= ggplot2::unit(0.5,'cm'),
    legend.background= ggplot2::element_blank(),
    #legend.box="bottom",
    legend.position ="bottom",legend.justification = c(1,0)
  )

#vol3
fc_clean =diff.result[is.finite(diff.result$log2FC), c("Label", "log2FC")]
fc_split = split(fc_clean$log2FC, fc_clean$Label)
x_range_per_label =lapply(fc_split, function(x) {
  x_abs <- abs(x)
  xmax <- max(x_abs, na.rm = TRUE) + 0.2
  xmin <- max(min(x_abs, na.rm = TRUE), 0.0001)
  c(xmin = xmin, xmax = xmax)
})
x_range_df = as.data.frame(do.call(rbind, x_range_per_label))
x_range_df$Label = rownames(x_range_df)
rownames(x_range_df) = NULL



get_curveData = function(xmin, xmax, fc_cutoff, p_cutoff) {
  xValues = seq(xmin, xmax, by = 0.0001)
  yValues = 1 / xValues + (-log10(p_cutoff))
  data.frame(
    xpos = c(xValues + log2(fc_cutoff), -(xValues +  log2(fc_cutoff))),
    ypos = c(yValues, yValues)
  )
}

curveData_list= lapply(1:nrow(x_range_df), function(i) {
  xmin_i =x_range_df$xmin[i]
  xmax_i = x_range_df$xmax[i]
  label_i=x_range_df$Label[i]
  df = get_curveData(xmin_i, xmax_i, cutoff_fc, cutoff_p)
  df$Label=label_i
  df
})

curveData_all =do.call(rbind, curveData_list)




plotData =diff.result  |> 
  dplyr::mutate(`-log10pvalue`=-log10(pvalue)) |> 
  dplyr::mutate(
    curveY = ifelse(log2FC > 0,
                    1 / (log2FC - log2(cutoff_fc)) + (-log10(cutoff_p)),
                    1 / (-log2FC - log2(cutoff_fc)) + (-log10(cutoff_p))),
    newDiff = dplyr::case_when(#Screening for more significant differences
      `-log10pvalue` > curveY & log2FC > log2(cutoff_fc) ~ 'up',
      `-log10pvalue` > curveY & log2FC < -log2(cutoff_fc) ~ 'down',
      TRUE ~ 'nosig'
    )
  ) |> 
  merge(x_range_df, by = "Label")







background.dat =data.frame(
  plotData |>  dplyr::group_by(Label) |> dplyr:: filter(log2FC>0) |>
    dplyr::summarise("y.localup"=max(log2FC)),
  plotData  |> dplyr:: group_by(Label) |> dplyr::filter(log2FC<0) |>
    dplyr::summarise("y.localdown"=min(log2FC)),
  x.local=seq(1:length(unique(plotData$Label)))
) |> dplyr:: select(-Label.1)


x.number =background.dat  |>   dplyr:: select(Label,x.local)
dat.plot = plotData |>   dplyr::left_join(x.number,by = "Label")

#selecting top-up and top-down proteins
dat.marked.up = plotData  |>  dplyr::filter(newDiff=="up")  |>
  dplyr::group_by(Label)  |> dplyr::arrange(-log2FC)  |>
  dplyr::top_n(top_marker,abs(log2FC))
dat.marked.down =plotData  |>  dplyr::filter(significance=="down") |>
  dplyr::group_by(Label)  |>  dplyr::arrange(log2FC)  |>
  dplyr:: top_n(top_marker,abs(log2FC))
dat.marked = dat.marked.up  |>  dplyr::bind_rows(dat.marked.down)
#referring group information data
dat.infor = background.dat  |>
  dplyr::mutate("y.infor"=rep(0,length(Label)))

down_up_counts=plotData |>
  dplyr::group_by(Label, newDiff) |>
  dplyr::summarise(cout = dplyr::n(), .groups = "drop") |>
  tidyr::complete(
    Label,
    newDiff = c("up", "down"),  
    fill = list(cout = 0)
  ) |> 
  dplyr:: filter(newDiff!="nosig")

dat.count= dat.infor|>
  dplyr::left_join(down_up_counts,by = "Label")



plt_vol3=
  ggplot2::ggplot()+
  ggplot2::geom_point(ggplot2::aes(x = log2FC, y = `-log10pvalue`),data = plotData  |>  dplyr::filter(newDiff == 'nosig') |>  dplyr::filter(is.finite(log2FC), is.finite(`-log10pvalue`)),size =2,color="grey70",alpha = 0.5)+
  ggplot2::geom_point(ggplot2::aes(x = log2FC, y = `-log10pvalue`,color=newDiff),data = plotData  |>  dplyr:: filter(newDiff!= 'nosig') |>  dplyr::filter(is.finite(log2FC), is.finite(`-log10pvalue`)),size =2,alpha = 0.6)+
  ggplot2::geom_line(data = curveData_all,  ggplot2::aes(x = xpos, y = ypos), lty = 2,  col = "grey10", lwd = 0.6) +
  ggrepel::geom_text_repel(data = plotData |> dplyr::arrange(dplyr::desc(abs(log2FC)))  |> 
                             dplyr::group_by(newDiff) |> 
                             dplyr::filter(!is.infinite(log2FC)) |> 
                             dplyr::filter(newDiff!="nosig") |> 
                             dplyr::slice_head(n = top_marker),
                           ggplot2::aes(x = log2FC, y = `-log10pvalue`, label = get(Feature)),
                           force = top_marker*1.5, color = 'black', size = 2,
                           point.padding = 0.5, hjust = 0.5,
                           arrow = ggplot2::arrow(length = grid::unit(0.02, "npc"),
                                         type = "open", ends = "last"),
                           segment.color="black",
                           segment.size = 0.3,
                           nudge_x = 0,
                           nudge_y = 1)+
  ggplot2:: geom_text(data = down_up_counts |> dplyr:: filter(newDiff == "up") ,
                      ggplot2::aes(x = Inf, y = Inf, label = paste0("n(up)=",cout)),
                      vjust = 2, hjust = 1, size = 3, color = "#B43F53", inherit.aes = FALSE) +
  ggplot2::geom_text(data = down_up_counts |> dplyr:: filter(newDiff == "down"),
                     ggplot2::aes(x = -Inf, y = Inf, label  = paste0("n(down)=",cout)),
                     vjust = 2, hjust = 0, size = 3, color = "#246367", inherit.aes = FALSE)+
  ggplot2:: labs(x = "log2(Fold change)",
                 y =paste0("-log10(pvalue)"),
                 title =paste0("Volcano of ",Feature) ) +
 
  ggplot2::facet_wrap(~Label,scales = "free")+
  ggplot2::coord_cartesian(ylim = c(0, max(plotData$`-log10pvalue`, na.rm = TRUE) * 1.5))+
  #ggplot2::ylim(0, max(plotData$`-log10pvalue`,na.rm = TRUE)*1.5)+
  ggplot2::scale_color_manual(values = c('up'="#B43F53", 
                                'down'='#246367'))+
  ggplot2::guides(color= ggplot2::guide_legend(title="Type"))+
 ggplot2::theme_bw(base_size = 10)+
  ggplot2::theme(plot.title = ggplot2::element_text(size = 15, hjust = 0.5),
                 plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5),
                 plot.background = ggplot2::element_rect(fill = "white", color = NA),
                 panel.background = ggplot2::element_rect(fill = "white", color = NA),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid = ggplot2::element_blank(),
                 panel.ontop = FALSE)



result=list(pltvol1=plt_vol,pltvol2=plt.vol2,plt.curve.vol=plt_vol3,
            dt_vol=diff.result |>as.data.frame() |>dplyr:: filter(significance!="nosig"),
            dt_vol_curve=plotData |> as.data.frame()|>dplyr:: filter(newDiff!="nosig") |> dplyr::rename(curve.significance=newDiff) |> dplyr::select(-xmin,-xmax)
            )

return(result)

}
