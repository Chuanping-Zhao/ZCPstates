#' @title uniprotEnrich_plot: Visualization for Enrichment Analysis
#' @description This function visualizes the results of enrichment analysis, including pathway or GO term enrichment, in bar plot, point plot, or network plot formats.
#'
#' @param dt A data frame containing enrichment analysis results. Must include columns such as `Pathway`, `p_value`, `FDR`, `counts`, `Enrich_factor`, etc.
#' @param plot.type A character string specifying the type of plot. Options are:
#' \itemize{
#'   \item `"bar"`: Bar plot for enrichment visualization.
#'   \item `"point"`: Point plot for enrichment visualization.
#'   \item `"network"`: Network plot showing relationships among enriched terms.
#' }
#' Default is `"point"`.
#' @param enrich.type A character string specifying the type of enrichment. Options are:
#' \itemize{
#'   \item `"pathway"`: For pathway enrichment analysis.
#'   \item `"GO"`: For Gene Ontology enrichment analysis.
#' }
#' Default is `"pathway"`.
#' @param GO.subset A character vector specifying the GO subset(s) to include in the plot. Options are:
#' \itemize{
#'   \item `"BP"`: Biological Process.
#'   \item `"CC"`: Cellular Component.
#'   \item `"MF"`: Molecular Function.
#'   \item `"All"`: All three GO subsets.
#' }
#' Default is `"CC"`.
#' @param network.fontsize Numeric value specifying font size for network labels. Default is `2`.
#' @param network.minclustersize Numeric value specifying the minimum size for clusters in network plots. Default is `2`.
#' @param network.simMethod A character string specifying the similarity method for clustering in network plots. Options are:
#' \itemize{
#'   \item `"jaccard"`: Jaccard similarity.
#'   \item `"cosine"`: Cosine similarity.
#'   \item `"cor"`: Correlation-based similarity.
#' }
#' Default is `"jaccard"`.
#'
#' @return A `ggplot` or `aPEAR` plot object, depending on the `plot.type` selected.
#'
#' @examples
#' # Example data
#'library(ZCPstates)
#' data("uniprotEnrichplot.demo.GO.CC", package = "ZCPstates")
#' go.cc.top10 <- as.data.frame(uniprotEnrichplot.demo.GO.CC) |> dplyr:: filter(Type=="CC") |>  dplyr:: slice_max(order_by = counts, n = 10)
#' # Plot pathway enrichment results as a bar plot
#' uniprotEnrich_plot(dt = go.cc.top10,plot.type = c("bar","point","network")[1], enrich.type=c("pathway","GO")[2],GO.subset=c("BP","CC","MF","All")[2])
#'
#'# Plot pathway enrichment results as a point plot
#'uniprotEnrich_plot(dt = go.cc.top10,plot.type = c("bar","point","network")[2], enrich.type=c("pathway","GO")[2],GO.subset=c("BP","CC","MF","All")[2])
#'
#' # Plot GO enrichment results as a network plot
#' uniprotEnrich_plot(dt = uniprotEnrichplot.demo.GO.CC, plot.type = c("bar","point","network")[3],enrich.type=c("pathway","GO")[2], GO.subset=c("BP","CC","MF","All")[2],network.fontsize=2,network.minclustersize=2,network.simMethod=c("jaccard", "cosine", "cor")[1])
#'
#' @importFrom dplyr select rename mutate filter arrange
#' @importFrom stringr str_extract str_remove_all str_remove str_trim
#' @importFrom ggplot2 ggplot aes geom_point geom_col geom_text scale_y_continuous coord_flip theme_bw theme labs scale_color_gradientn scale_fill_gradientn facet_wrap
#' @importFrom crayon yellow
#' @export
#'
uniprotEnrich_plot=function(
    dt,
    plot.type=c("bar","point","network")[2],
    enrich.type=c("pathway","GO")[1],
    GO.subset=c("BP","CC","MF","All")[2],
    network.fontsize=2,
    network.minclustersize=2,
    network.simMethod=c("jaccard", "cosine", "cor")[1]
){

  switch (enrich.type,
          "pathway" = {
            dt.plot=as.data.frame(dt) |>
              dplyr::select(any_of(c("Type","Pathway","p_value","FDR","hited.proteins","counts","Enrich_factor"))) |>
              dplyr::rename("ONTOLOGY"="Type",
                            "Description"="Pathway",
                            "pvalue"="p_value",
                            "p.adjust"="FDR",
                            "geneID"="hited.proteins",
                            "Count"="counts",
                            "NES"="Enrich_factor"
              ) |>
              dplyr::mutate(ID=Description)
          },
          "GO"={
            if(!any(GO.subset%in%c("BP","CC","MF"))){
              GO.subset <- c("BP","CC","MF")
              cat(crayon::yellow("Defaulting to 'All' GO enrichment data.--zcp\n"))
            }

            if(length(GO.subset)==1){
              GO.subset <- GO.subset[[1]]
              dt.plot=as.data.frame(dt) |> dplyr:: mutate(
                ID = stringr::str_extract(Pathway, "\\[.*?\\]") |> stringr::str_remove_all("\\[|\\]"),
                Description = stringr::str_remove(Pathway, "\\[.*?\\]") |> stringr::str_trim()
              ) |>
                dplyr::select(!Pathway) |>
                dplyr::rename("ONTOLOGY"="Type",
                              "NES"="Enrich_factor",
                              "pvalue"="p_value",
                              "p.adjust"="FDR",
                              "geneID"="hited.proteins",
                              "Count"="counts"
                )   |>
                dplyr::filter(ONTOLOGY==GO.subset)
            }else{
              dt.plot=as.data.frame(dt) |> dplyr:: mutate(
                ID = stringr::str_extract(Pathway, "\\[.*?\\]") |> stringr::str_remove_all("\\[|\\]"),
                Description = stringr::str_remove(Pathway, "\\[.*?\\]") |> stringr::str_trim()
              ) |>
                dplyr::select(!Pathway) |>
                dplyr::rename("ONTOLOGY"="Type",
                              "NES"="Enrich_factor",
                              "pvalue"="p_value",
                              "p.adjust"="FDR",
                              "geneID"="hited.proteins",
                              "Count"="counts"
                )
            }


          }
  )






  switch (plot.type,
          "network" = {#
            plt.network=suppressMessages(
              suppressWarnings(
                enrichmentNetwork(dt.plot,
                                         drawEllipses = TRUE,
                                         simMethod=network.simMethod,#方法参数

                                         repelLabels=TRUE,
                                         colorBy = 'NES',
                                         nodeSize = 'Count',
                                         colorType = c("nes", "pval")[1],
                                         minClusterSize=network.minclustersize,
                                         fontSize = network.fontsize,
                                         verbose = FALSE)+
                  #ggplot2::scale_color_gradientn(colours = c("#006AD2","white","#AF217C"),name = "Enrich factor")+
                  ggplot2::scale_color_gradientn(colours = c("#1A5592","white","#B83D3D"),name = "NES")+
                  #viridis::scale_color_viridis(option = viridis.color,name="Enrich factor")+#,direction = -1
                  ggplot2::guides(size = ggplot2::guide_legend(title = "Pathway size"))+
                  ggplot2:: theme(legend.text = ggplot2::element_text(size = 10),
                                  legend.position = "left")
              )
              )
          },
          "point"={#
            plt.network=
              ggplot2::ggplot(data = dt.plot,ggplot2::aes(x=Description,y=Count,size=Count,color=-log10(pvalue)))+
              ggplot2::geom_point(alpha=0.9)+
              ggplot2::scale_size_continuous(range = c(2, 5), name = "Count") +
              ggplot2::facet_wrap(~ONTOLOGY,scales = "free")+
              ggplot2::theme_bw(base_size = 8)+
              ggplot2::theme(
                panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank()#,
                # legend.position = "bottom"
              )+
              #ggplot2::labs(x="Count")+
              ggplot2::coord_flip()+
              #viridis::scale_color_viridis(option = "A",name="-log10 pvalue")+
              ggplot2::scale_color_gradientn(colours = c("#1A5592","white","#B83D3D"),name = "-log10 pvalue")+
              ggplot2::guides(size = "none")

          },
          "bar"={#
            dt.plot1=dt.plot |>  dplyr::arrange(Count) |>  dplyr::mutate(Description=factor(Description,Description)) |>
              dplyr::mutate(AdjustedCount = Count * 0.8)

            plt.network=
              ggplot2::ggplot(data = dt.plot1,ggplot2::aes(x=Description,y=Count,fill=-log10(pvalue)))+
              ggplot2::geom_col(color="black",linewidth=0.3)+
              ggplot2::geom_text(
                ggplot2::aes(label = paste0("NES: ", sprintf("%0.2f", NES)), y = Count),
                parse = TRUE,
                hjust = -0.1,
                size = 2
              )+
              ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, max(dt.plot1$Count) * 1.2))+
              ggplot2::facet_wrap(~ONTOLOGY,scales = "free")+
              ggplot2::theme_bw(base_size = 8)+
              ggplot2::theme(
                panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank()#,
                # legend.position = "bottom"
              )+

              ggplot2::coord_flip()+
              ggplot2::labs(y="Count",y=NULL)+
              ggplot2::scale_fill_gradientn(colours = c("#1A5592","white","#B83D3D"),name = "-log10 pvalue")

          }
  )
return(plt.network)
}
