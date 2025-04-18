% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jaccard_indices.R
\name{jaccard_indices}
\alias{jaccard_indices}
\title{Compute and Visualize Pairwise Jaccard Indices Between Sample Groups}
\usage{
jaccard_indices(
  dt,
  intensity.col = "Abundance",
  protein.col = "Protein",
  group.col = "BioReplicate",
  fillcolor_high = "#EC7607",
  fillcolor_low = "#077DEC",
  fontsize = 2
)
}
\arguments{
\item{dt}{A data.frame or tibble containing at least protein identifiers, group labels, and intensity values.}

\item{intensity.col}{Character. Column name that stores intensity or abundance values. Default is \code{"Abundance"}.}

\item{protein.col}{Character. Column name that stores protein identifiers. Default is \code{"Protein"}.}

\item{group.col}{Character. Column name that stores sample group or replicate names. Default is \code{"BioReplicate"}.}

\item{fillcolor_high}{Character. Color used for high Jaccard values in heatmap. Default is \code{"#EC7607"}.}

\item{fillcolor_low}{Character. Color used for low Jaccard values in heatmap. Default is \code{"#077DEC"}.}

\item{fontsize}{Numeric. Font size for the Jaccard value labels shown on the heatmap tiles. Default is \code{2}.}
}
\value{
A named list with two elements:
\describe{
\item{\code{plot}}{A ggplot object showing the lower triangle of the Jaccard index matrix.}
\item{\code{jaccard_indices}}{A data.frame in long format representing the lower triangle of the Jaccard matrix (columns: Var1, Var2, value).}
}
}
\description{
This function computes the pairwise Jaccard similarity indices of protein presence
across different groups (e.g., biological replicates), and visualizes the lower triangle
of the resulting Jaccard matrix as a heatmap.
}
\details{
It removes entries from the input data where the intensity column is missing, zero,
"0", empty string, or NaN before computation.
}
\examples{
dt <- data.frame(
  Protein = paste0("P", 1:10),
  BioReplicate = rep(c("Rep1", "Rep2"), each = 5),
  Abundance = c(runif(5, 5, 10), runif(5, 0, 10))
)
res <- jaccard_indices(dt,
intensity.col = "Abundance",
protein.col="Protein",
group.col="BioReplicate",
fillcolor_high="#EC7607",
fillcolor_low="#077DEC",
fontsize=2)
res$plot  #View the heatmap

}
