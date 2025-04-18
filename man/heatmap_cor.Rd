% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmap_cor.R
\name{heatmap_cor}
\alias{heatmap_cor}
\title{Plot Correlation Matrix Heatmap with Annotations}
\usage{
heatmap_cor(
  x,
  cor.method = "spearman",
  Type = "Full",
  highcolor = "#492952",
  lowcolor = "#1e4668",
  midpoints = 0.5,
  lowpoints = 0,
  highpoints = 1,
  Title = NULL,
  geom_tile_color = "white",
  geom_tile_text_size = 0.3,
  fontcolor = "white",
  lab.x.fontsize = 6,
  lab.y.fontsize = 6,
  cor.fontsize = 2
)
}
\arguments{
\item{x}{A numeric matrix or data frame (e.g., protein abundances with samples in columns).}

\item{cor.method}{Character. Correlation method to use. Options: \code{"pearson"}, \code{"kendall"}, or \code{"spearman"}. Default is \code{"spearman"}.}

\item{Type}{Character. Type of matrix to display: \code{"Full"}, \code{"upper_tri"}, or \code{"lower_tri"}. Default is \code{"Full"}.}

\item{highcolor}{Character. Color for high correlation values (close to 1). Default is \code{"#492952"}.}

\item{lowcolor}{Character. Color for low correlation values (close to 0). Default is \code{"#1e4668"}.}

\item{midpoints}{Numeric. Midpoint value for color scale. Default is \code{0.5}.}

\item{lowpoints}{Numeric. Lower bound of color scale. Default is \code{0}.}

\item{highpoints}{Numeric. Upper bound of color scale. Default is \code{1}.}

\item{Title}{Character. Plot title. Default is \code{NULL}.}

\item{geom_tile_color}{Character. Color of the tile border. Default is \code{"white"}.}

\item{geom_tile_text_size}{Numeric. Border line size of tiles. Default is \code{0.3}.}

\item{fontcolor}{Character. Color of text labels on heatmap. Default is \code{"white"}.}

\item{lab.x.fontsize}{Numeric. Font size for x-axis labels. Default is \code{6}.}

\item{lab.y.fontsize}{Numeric. Font size for y-axis labels. Default is \code{6}.}

\item{cor.fontsize}{Numeric. Font size for correlation value labels. Default is \code{2}.}
}
\value{
A named list with:
\describe{
\item{\code{plot}}{A \code{ggplot2} object showing the correlation heatmap.}
\item{\code{correlation}}{The numeric correlation matrix used for plotting.}
}
}
\description{
This function computes the correlation matrix between numeric columns of an input matrix or data frame,
and plots a customizable heatmap with correlation values labeled on tiles. It supports showing the full,
upper, or lower triangle of the correlation matrix and allows user-defined color schemes, font sizes,
and placement of legends.
}
\examples{
# Example
mat <- matrix(rnorm(100), ncol = 10)
colnames(mat) <- paste0("Sample", 1:10)
res <- heatmap_cor(mat, cor.method = "pearson", Type = "lower_tri")
res$plot  # Show the plot

}
