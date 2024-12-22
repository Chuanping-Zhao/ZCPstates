#' @title uniprotEnrich.demo.library.uniprot: Example UniProt Database
#' @description This dataset contains a sample UniProt database for proteome UP000006548.
#' It includes columns such as "Entry", "Pathway", and Gene Ontology annotations.
#' @format A data frame with N rows and M columns.
#' \describe{
#'   \item{Entry}{Protein ID from UniProt.}
#'   \item{Pathway}{Pathway annotations associated with proteins.}
#'   \item{...}{Additional columns for GO annotations.}
#' }
#' @source \url{https://uniprot.org}
"uniprotEnrich.demo.library.uniprot"

#' @title uniprotEnrich.demo.diff: Example Dataset for Enrichment analysis
#' @description This dataset contains example data that has been finalized after differential analysis and is intended for downstream enrichment analysis.
#' It includes protein IDs and differential conditions.
#' @format  A data frame with N rows and M columns.
#' \describe{
#'   \item{ProteinIDs}{Protein identifiers.}
#'   \item{sig}{Differential condition marker (e.g., "Up", "Down").}
#'   \item{...}{Additional columns as required.}
#' }
#' @source Example data
"uniprotEnrich.demo.diff"

#'@title uniprotEnrichplot.demo.GO.CC: Example Dataset for visualing the enrichment analysis result
#'@description This dataset is derived from the analysis results of the uniprotEnrich function, specifically focusing on the CC (cellular component) subset of GO enrichment. It is intended to demonstrate the visualization capabilities of the uniprotEnrich_plot function.
#'@format  A data frame with N rows and M columns.
#'#' \describe{
#'   \item{Type}{GO subterm identifiers.}
#'   \item{counts}{Enriched pathway size.}
#'    \item{Enrich_factor}{NES value for each pathway.}
#'   \item{...}{Additional columns as required.}
#' }
#'@source Example data
"uniprotEnrichplot.demo.GO.CC"
