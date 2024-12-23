% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uniprotEnrich_plot.R
\name{uniprotEnrich_plot}
\alias{uniprotEnrich_plot}
\title{uniprotEnrich_plot: Visualization for Enrichment Analysis}
\usage{
uniprotEnrich_plot(
  dt,
  plot.type = c("bar", "point", "network")[2],
  enrich.type = c("pathway", "GO")[1],
  GO.subset = c("BP", "CC", "MF", "All")[2],
  network.fontsize = 2,
  network.minclustersize = 2,
  network.simMethod = c("jaccard", "cosine", "cor")[1]
)
}
\arguments{
\item{dt}{A data frame containing enrichment analysis results. Must include columns such as \code{Pathway}, \code{p_value}, \code{FDR}, \code{counts}, \code{Enrich_factor}, etc.}

\item{plot.type}{A character string specifying the type of plot. Options are:
\itemize{
\item \code{"bar"}: Bar plot for enrichment visualization.
\item \code{"point"}: Point plot for enrichment visualization.
\item \code{"network"}: Network plot showing relationships among enriched terms.
}
Default is \code{"point"}.}

\item{enrich.type}{A character string specifying the type of enrichment. Options are:
\itemize{
\item \code{"pathway"}: For pathway enrichment analysis.
\item \code{"GO"}: For Gene Ontology enrichment analysis.
}
Default is \code{"pathway"}.}

\item{GO.subset}{A character vector specifying the GO subset(s) to include in the plot. Options are:
\itemize{
\item \code{"BP"}: Biological Process.
\item \code{"CC"}: Cellular Component.
\item \code{"MF"}: Molecular Function.
\item \code{"All"}: All three GO subsets.
}
Default is \code{"CC"}.}

\item{network.fontsize}{Numeric value specifying font size for network labels. Default is \code{2}.}

\item{network.minclustersize}{Numeric value specifying the minimum size for clusters in network plots. Default is \code{2}.}

\item{network.simMethod}{A character string specifying the similarity method for clustering in network plots. Options are:
\itemize{
\item \code{"jaccard"}: Jaccard similarity.
\item \code{"cosine"}: Cosine similarity.
\item \code{"cor"}: Correlation-based similarity.
}
Default is \code{"jaccard"}.}
}
\value{
A \code{ggplot} or \code{aPEAR} plot object, depending on the \code{plot.type} selected.
}
\description{
This function visualizes the results of enrichment analysis, including pathway or GO term enrichment, in bar plot, point plot, or network plot formats.
}
\examples{
# Example data
library(ZCPstates)
data("uniprotEnrichplot.demo.GO.CC", package = "ZCPstates")
go.cc.top10 <- as.data.frame(uniprotEnrichplot.demo.GO.CC) |> dplyr:: filter(Type=="CC") |>  dplyr:: slice_max(order_by = counts, n = 10)
# Plot pathway enrichment results as a bar plot
uniprotEnrich_plot(dt = go.cc.top10,plot.type = c("bar","point","network")[1], enrich.type=c("pathway","GO")[2],GO.subset=c("BP","CC","MF","All")[2])

# Plot pathway enrichment results as a point plot
uniprotEnrich_plot(dt = go.cc.top10,plot.type = c("bar","point","network")[2], enrich.type=c("pathway","GO")[2],GO.subset=c("BP","CC","MF","All")[2])

# Plot GO enrichment results as a network plot
uniprotEnrich_plot(dt = uniprotEnrichplot.demo.GO.CC, plot.type = c("bar","point","network")[3],enrich.type=c("pathway","GO")[2], GO.subset=c("BP","CC","MF","All")[2],network.fontsize=2,network.minclustersize=2,network.simMethod=c("jaccard", "cosine", "cor")[1])

}
