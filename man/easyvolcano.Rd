% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/easyvolcano.R
\name{easyvolcano}
\alias{easyvolcano}
\title{Visualize easyDiff Results with Multiple Volcano Plots}
\usage{
easyvolcano(
  diff,
  cutoff_fc = 2,
  cutoff_p = 0.05,
  top_marker = 0,
  max_overlaps = 0,
  Feature = "Protein",
  use_adj_p = FALSE
)
}
\arguments{
\item{diff}{The result of an \code{easydiff} analysis, provided as a tibble or data frame with columns:
\code{Label}, \code{Protein}, \code{log2FC}, and \code{pvalue}. Each row represents a feature's result from one comparison.}

\item{cutoff_fc}{A numeric threshold for the fold-change cutoff (on the log2 scale). Default is \code{2}.}

\item{cutoff_p}{A numeric threshold for the p-value significance cutoff. Default is \code{0.05}.}

\item{top_marker}{Integer, number of top up- and down-regulated features to label per group. Default is \code{0}.}

\item{max_overlaps}{Maximum number of overlapping labels allowed, passed to \code{ggrepel::geom_text_repel()}. Default is \code{0}.}

\item{Feature}{Character string to describe the feature type (e.g., \code{"protein"}, \code{"metabolite"}) used in plot annotations. Default is \code{"Protein"}.}

\item{use_adj_p}{Logical. If \code{TRUE}, use the column \code{adj.pvalue} instead of \code{pvalue} for determining statistical significance. Default is \code{FALSE}.}
}
\value{
A list containing:
\describe{
\item{\code{pltvol1}}{Standard volcano plot with color-coded points by significance and group-wise annotations.}
\item{\code{pltvol2}}{Categorical bar-style volcano plot, summarizing significance across groups with optional top labels.}
\item{\code{plt.curve.vol}}{Curve-based volcano plot with dynamic significance curves per group and labeled top features.}
\item{\code{dt_vol}}{Filtered data frame of significant features based on fold change and (adjusted) p-value cutoffs.}
\item{\code{dt_vol_curve}}{Filtered data frame based on dynamic significance curves, including assigned \code{curve.significance}.}
}
}
\description{
\code{easyvolcano} generates volcano plots to visualize differential analysis results from \code{easydiff}.
It highlights features (e.g., proteins, peptides, metabolites) that show significant changes between conditions,
based on specified fold-change and significance thresholds, and offers three distinct visual styles for enhanced interpretation.
}
\examples{
\dontrun{
# Generate volcano plots with default settings:
volcano_output <- easyvolcano(diff, cutoff_fc = 2, cutoff_p = 0.05, top_marker = 10, Feature = "protein")

# Display each plot:
volcano_output$pltvol1
volcano_output$pltvol2
volcano_output$plt.curve.vol

# Examine significant results:
head(volcano_output$dt_vol)
head(volcano_output$dt_vol_curve)
}

}
