% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/easystringEnrich.R
\name{easystringEnrich}
\alias{easystringEnrich}
\title{easystringEnrich: Protein Enrichment Analysis Function}
\usage{
easystringEnrich(
  dt,
  proteincol.id = "Accessions",
  sep.pattern = ";",
  savepath = "enrichoutputfile",
  backgrounds = NULL,
  cutoff = 2,
  protein_aliases_path,
  protein_Info_path,
  protein_enrichment_terms_path
)
}
\arguments{
\item{dt}{A data table (or data frame) containing protein ID information.}

\item{proteincol.id}{A string specifying the column name that stores protein IDs. Default is "Accessions".}

\item{sep.pattern}{A string specifying the delimiter used to separate protein IDs. Default is ";".}

\item{savepath}{A string specifying the directory where the results and plots will be saved. Default is "enrichoutputfile".}

\item{backgrounds}{Optional background protein data table; if NULL, all proteins in the background will be used.}

\item{cutoff}{Numeric value; proteins with counts below this threshold will be filtered out. Default is 2.}

\item{protein_aliases_path}{A string providing the local file path for the protein aliases file
(e.g., "9606.protein.aliases.v12.0.txt.gz").}

\item{protein_Info_path}{A string providing the local file path for the protein info file
(e.g., "9606.protein.info.v12.0.txt.gz").}

\item{protein_enrichment_terms_path}{A string providing the local file path for the protein enrichment terms file
(e.g., "9606.protein.enrichment.terms.v12.0.txt.gz").}
}
\value{
A list containing:
\describe{
\item{enrichment}{The enrichment analysis results as a tibble.}
\item{gene_pathway_matrix}{The gene-pathway mapping matrix as a tibble.}
\item{enrichplots}{A named list of ggplot objects, each corresponding to a pathway category.}
}
}
\description{
Performs enrichment analysis on the input protein data. This function integrates protein aliases, protein information,
and protein enrichment terms data (string database v12.0), calculates p-values and FDR using the hypergeometric test, and generates pathway
enrichment plots.
}
