#' Plot an UpSet Plot for Protein Group Memberships
#'
#' This function generates an UpSet plot using ggupset and ggplot2 to visualize the distribution of proteins
#' (e.g., peptide sequences) across different sample groups or conditions.
#'
#' @param longdata A data.frame or data.table in long format, typically with at least one column for protein IDs
#'   and one column representing group membership (e.g., sample or condition names).
#' @param group.by A string specifying the column name to use for group membership (e.g., `"R.FileName"`).
#' @param proteincol.id A string specifying the column name that contains the unique protein or peptide identifiers
#'   (e.g., `"PEP.StrippedSequence"`).
#' @param top_n Integer. Number of top combinations (intersections) to display in the plot. Default is 10.
#' @param combmatrix_point_color Color for the combination matrix points. Default is "#0863A7".
#' @param combmatrix_line_color Color for the combination matrix lines. Default is "#8FD3BD".
#' @param combmatrix_linesize Numeric. Size (thickness) of the combination matrix lines. Default is 2.
#' @param alpha Numeric. Transparency for the bars. Range from 0 to 1. Default is 1 (fully opaque).
#'
#' @return A ggplot2 object representing the UpSet plot.
#'
#' @importFrom data.table as.data.table setnames
#' @importFrom ggplot2 ggplot aes_string geom_bar geom_text theme_classic after_stat
#' @export
#'
#' @examples
#' # Create demo data
#' set.seed(123)
#' demo_data <- data.frame(
#'   R.FileName = sample(paste0("Sample", 1:5), 100, replace = TRUE),
#'   PEP.StrippedSequence = sample(paste0("PEP", 1:20), 100, replace = TRUE),
#'   PG.ProteinAccessions = sample(c("P001", "P002", "P003"), 100, replace = TRUE),
#'   PG.ProteinDescriptions = sample(c("Kinase", "Transporter", "Ligase"), 100, replace = TRUE),
#'   PG.ProteinNames = sample(c("KIN1_HUMAN", "TRP1_HUMAN", "LIG1_HUMAN"), 100, replace = TRUE),
#'   PG.Coverage = paste0(sample(5:30, 100, replace = TRUE), "%"),
#'   PG.IsSingleHit = sample(c(TRUE, FALSE), 100, replace = TRUE),
#'   PG.Qvalue = runif(100),
#'   PG.Quantity = runif(100, 1e3, 1e5)
#' )
#'
#' # Plot using easy_upset
#' easy_upset(
#'   longdata = demo_data,
#'   group.by = "R.FileName",
#'   proteincol.id = "PEP.StrippedSequence",
#'   top_n = 8,
#'   combmatrix_point_color = "#4C72B0",
#'   combmatrix_line_color = "#55A868",
#'   alpha = 0.8
#' )
#'
easy_upset = function(
    longdata,
    group.by,
    proteincol.id,
    top_n = 10,
    combmatrix_point_color= "#0863A7",
    combmatrix_line_color=  "#8FD3BD",
    combmatrix_linesize=2,
    alpha = 1
) {

  viridis_colors = viridis::viridis(top_n)
  dt= data.table::as.data.table(longdata)
  dt.upset = dt[, .(group_list = list(get(group.by))), by = proteincol.id]
  data.table::setnames(dt.upset, "group_list", group.by)


  p = ggplot2::ggplot(data = base::as.data.frame(dt.upset),
                       ggplot2::aes_string(x = group.by)) +
    ggupset::scale_x_upset(order_by = "freq", n_intersections = top_n) +
    ggplot2::geom_bar(fill = viridis_colors, color = "black", alpha = alpha) +
    ggplot2::geom_text(stat = "count", ggplot2::aes(label = ggplot2::after_stat(count)), vjust = -0.5) +
    ggplot2::theme_classic() +
    ggupset::theme_combmatrix(
      combmatrix.panel.point.color.fill = combmatrix_point_color,
      combmatrix.panel.line.color = combmatrix_line_color,
      combmatrix.panel.line.size = combmatrix_linesize
    )

  return(p)
}
