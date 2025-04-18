% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boxdensity.R
\name{boxdensity}
\alias{boxdensity}
\title{Plot Boxplot and Density Plot for Quantitative QC}
\usage{
boxdensity(
  dt,
  intensity.col = "Abundance",
  sample.col = "BioReplicate",
  box.outlier.color = "#70AD47",
  box.fill = "#00B050",
  density.fill = "#00B050",
  logTrans = TRUE
)
}
\arguments{
\item{dt}{A data.frame or tibble containing the intensity values and sample information.}

\item{intensity.col}{Character. Name of the column representing intensity values. Default is \code{"Abundance"}.}

\item{sample.col}{Character. Name of the column representing sample or replicate labels. Default is \code{"BioReplicate"}.}

\item{box.outlier.color}{Character. Color of outliers in boxplot. Default is \code{"#70AD47"}.}

\item{box.fill}{Character. Fill color of the boxplot. Default is \code{"#00B050"}.}

\item{density.fill}{Character. Fill color of the density plot. Default is \code{"#00B050"}.}

\item{logTrans}{Logical. Whether to apply log2 transformation to the intensity column. If \code{FALSE}, log2 will be applied. Default is \code{TRUE}.}
}
\value{
A \code{ggplot} patchwork object combining a density plot and a boxplot.
}
\description{
This function creates a vertically stacked plot composed of a density plot and a boxplot
to visualize distribution and variability of intensity (e.g. protein abundance) across replicates.
}
\examples{
# Example data
dt <- data.frame(
  Protein = paste0("P", 1:20),
  BioReplicate = rep("QC1", 20),
  Abundance = runif(20, 1, 1000)
)
boxdensity(dt, logTrans = FALSE)

}
