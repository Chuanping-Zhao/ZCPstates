% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_comparison_matrix.R
\name{generate_comparison_matrix}
\alias{generate_comparison_matrix}
\title{Generate All Pairwise Comparison Matrix}
\usage{
generate_comparison_matrix(
  sampleInfo,
  group_col = "Condition",
  force_denominator = NULL
)
}
\arguments{
\item{sampleInfo}{A data frame or tibble containing at least one column specifying group information.}

\item{group_col}{Character. The name of the column in \code{sampleInfo} that contains the condition or group labels. Default is \code{"Condition"}.}

\item{force_denominator}{Optional character vector. Conditions specified here will always be placed in the denominator (i.e., assigned \code{-1} in the comparison vector) when involved in a comparison.}
}
\value{
A matrix with rows representing comparisons (e.g., \code{"A_vs_B"}) and columns corresponding to each unique condition. Each row is a contrast vector with \code{1} for numerator, \code{-1} for denominator, and \code{0} otherwise.
}
\description{
This function creates a comparison matrix for all pairwise combinations of experimental conditions,
useful for differential analysis. It supports forcing specific groups to always appear in the denominator.
}
\examples{
sampleInfo <- tibble::tibble(
  Run = paste0("Run", 1:6),
  Condition = c("A", "A", "B", "B", "C", "C")
)

generate_comparison_matrix(sampleInfo)
generate_comparison_matrix(sampleInfo, force_denominator = c("C"))

}
