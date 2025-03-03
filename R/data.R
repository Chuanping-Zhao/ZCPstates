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

#' @title Glycan_network_demo: Example Glycan Network Data
#'
#' @description
#' This dataset provides a sample glycan network for demonstration, containing
#' columns such as \code{Annotated Sequence}, \code{Glycan composition}, and
#' \code{Master Protein Accessions}. It can be used to test or illustrate
#' glycan-related functions in the package.
#'
#' @format A data frame with 21 rows and 3 columns:
#' \describe{
#'   \item{Annotated Sequence}{The annotated peptide sequence.}
#'   \item{Glycan composition}{The glycan composition, e.g. \code{HexNAc(2)Hex(8)}.}
#'   \item{Master Protein Accessions}{Protein accession identifiers (e.g., \code{Q8TCJ2}).}
#' }
#'
#' @details
#' Each row corresponds to a particular glycopeptide with its associated glycan
#' composition and the protein accession from which it originates. This dataset
#' can be loaded with:
#' \preformatted{
#'   data("Glycan_network_demo", package = "ZCPstates")
#' }
#'@source Example data
"Glycan_network_demo"
