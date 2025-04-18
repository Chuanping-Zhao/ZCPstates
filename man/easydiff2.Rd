% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/easydiff2.R
\name{easydiff2}
\alias{easydiff2}
\title{Perform differential expression analysis on a wide-format protein matrix}
\usage{
easydiff2(
  data = protein_wide,
  contrast.matrix = comparison,
  groupInfo = sampleInfo,
  NAfilter.cutoff = NULL,
  save_fitted_models = TRUE,
  protein.col = "x",
  log2Trans = FALSE
)
}
\arguments{
\item{data}{A wide-format data.frame or tibble of protein expression values. Rows are proteins, columns are sample Runs. One column must contain protein IDs.}

\item{contrast.matrix}{A contrast matrix specifying comparisons to make between experimental groups. Columns must match group names in \code{groupInfo$Condition}.}

\item{groupInfo}{A data.frame or tibble with columns \code{Run}, \code{Condition}, and \code{BioReplicate}, specifying the sample group assignments.}

\item{NAfilter.cutoff}{Optional. If specified, proteins with missing percentage above this threshold will be filtered out. (e.g. 0.5)}

\item{save_fitted_models}{Logical. If \code{TRUE}, optionally retain fitted model objects per protein. Default is \code{TRUE}.}

\item{protein.col}{Column name containing protein IDs in the input \code{data}. Default is \code{"x"}.}

\item{log2Trans}{Logical. If \code{FALSE}, input data will be log2-transformed. Default is \code{FALSE}.}
}
\value{
A data.table of differential expression results across all proteins, including:
\itemize{
\item \code{Protein}: Protein ID
\item \code{Label}: Contrast label (e.g. "A_vs_B")
\item \code{log2FC}: Estimated log2 fold change from linear model
\item \code{SE}: Standard error
\item \code{Tvalue}: t statistic
\item \code{DF}: Residual degrees of freedom
\item \code{pvalue}: Raw p-value
\item \code{adj.pvalue}: Multiple testing adjusted p-value (BH method)
}
}
\description{
This function performs linear model-based group comparison on log2-intensity
proteomics data. It automatically reshapes the input data, merges sample annotations,
handles missing value filtering, and fits fixed-effect linear models for each protein.
}
\details{
The function uses \code{lm(ABUNDANCE ~ GROUP)} as the base model. Missing values are optionally filtered. Contrasts are specified manually.
Internally reshapes the data, merges group info, and splits by protein for independent model fitting.
}
\examples{
# example
library(ZCPstates)
protein_wide <- data.frame(
x = c("P001", "P002", "P003"),
Sample1 = c(12.3, 13.1, 14.0),
Sample2 = c(12.5, 13.0, 14.2),
Sample3 = c(11.8, 12.9, 13.9),
Sample4 = c(14.1, 13.5, 15.0),
Sample5 = c(13.9, 13.7, 14.9),
Sample6 = c(14.2, 13.6, 15.1)
)
sampleInfo <- data.frame(
Run = paste0("Sample", 1:6),
Condition = c("A", "A", "A", "B", "B", "B"),
BioReplicate = paste0("R", 1:6)
)
comparison=generate_comparison_matrix(sampleInfo,
group_col = "Condition",
force_denominator = NULL)

result <- easydiff2(
  data = protein_wide,
  contrast.matrix = comparison,
  groupInfo = sampleInfo,
  NAfilter.cutoff = 0.5,
  save_fitted_models = TRUE,
  protein.col = "x",
  log2Trans = TRUE
)
head(result)
}
