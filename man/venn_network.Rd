% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/venn_network.R
\name{venn_network}
\alias{venn_network}
\title{Plot a Set-Based Venn Network with Side Tables}
\usage{
venn_network(
  df,
  protein_col = "Protein",
  method_col = "Method",
  max_show = 10,
  font_size = 10,
  font_face = "plain",
  font_family = "mono",
  alpha = 0.3,
  show_labels = FALSE,
  show_table = FALSE,
  set_colors = NULL
)
}
\arguments{
\item{df}{A data.frame or data.table.}

\item{protein_col}{Column name for protein IDs.}

\item{method_col}{Column name for method/group.}

\item{max_show}{Max proteins per side table.}

\item{font_size}{Font size for side tables.}

\item{font_face}{Font face (e.g., "plain", "bold").}

\item{font_family}{Font family for all text.}

\item{alpha}{Table background alpha (0–1).}

\item{show_labels}{Logical. Whether to show protein node labels.}

\item{show_table}{Logical. Whether to show unique element side tables.}

\item{set_colors}{Optional named vector of colors for each set (names must match methods).}
}
\value{
A patchwork plot with network and side tables (if \code{show_table = TRUE}), or network only.
}
\description{
This function visualizes overlapping sets as a network graph, with individual
elements and their memberships. Unique elements are also shown in side tables.
}
\examples{
set.seed(123)
demo_proteins <- paste0("P", sprintf("\%05d", sample(10000:99999, 600)))
group1 <- sample(demo_proteins, 300)
group2 <- sample(demo_proteins, 300)
df_demo <- data.frame(
  Protein = c(group1, group2),
  Method = c(rep("A", length(group1)), rep("B", length(group2)))
)

venn_network(df_demo, show_labels = TRUE, show_table = TRUE)

venn_network(df_demo, show_table = TRUE, set_colors = c(A = "#3A7A77", B = "#60C080"))

}
