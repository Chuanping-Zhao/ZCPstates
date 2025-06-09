#' Plot Step Line Chart
#'
#' This function draws a step line chart with midpoint steps and value labels.
#' It supports grouping by method (e.g., "A", "B"), legend customization, and color/line width control.
#'
#' @param dt A data frame in long format with three columns: `Sample`, `Group`, and `count`.
#'   - `Sample`: Sample names (e.g., "s1", "s2", ...)
#'   - `Group`: Grouping factor (e.g., "A", "B")
#'   - `count`: Numeric values to plot
#' @param color_map A named character vector specifying colors for each group. Default is c("A" = "#C0504D", "B" = "#4F81BD").
#' @param line_width Numeric value specifying the width of the step lines. Default is 0.8.
#' @param x.lab Character of xlabs.
#'  @param y.lab Character of ylabs.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' dt <- data.frame(
#'   Sample = rep(paste0("s", 1:7), each = 2),
#'   Group = rep(c("A", "B"), times = 7),
#'    count = c(25,100, 31, 108,66,115, 10,80, 28,90, 40,94, 50,130 )
#' )
#' easy_stepline(dt)
#' easy_stepline(dt, color_map = c("darkred", "steelblue"), line_width = 1.2)
#'
#' @export
easy_stepline <- function(dt,
                          color_map = c("#C0504D",  "#4F81BD"),
                          x.lab="sample",
                          y.lab="counts",
                          line_width = 0.8) {

  df <- dt |>
    dplyr::mutate(Sample = factor(Sample, levels = unique(Sample))) |>
    dplyr::mutate(x = as.numeric(Sample) - 1)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = count, color = Group)) +
    ggplot2::geom_step(direction = "mid", linewidth = line_width) +
    ggplot2::geom_text(
      ggplot2::aes(label = count, y = count),
      size = 2.8, show.legend = FALSE, vjust = -0.5, hjust = 0.5
    ) +
    ggplot2::scale_x_continuous(
      breaks = unique(df$x),
      labels = unique(df$Sample),
      expand = c(0, 0.4)
    ) +
    ggplot2::scale_color_manual(values = color_map) +
    ggplot2::labs(x = x.lab, y = y.lab) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 15, hjust = 0.5),
      panel.grid = ggplot2::element_blank(),
      legend.position = c(0, 1),
      legend.justification = c(0, 1),
      legend.background = ggplot2::element_rect(
        fill = "#FFFFFFAA",
        colour = "black",
        linewidth = 0.3,
        linetype = "solid"
      ),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9)
    )

  return(p)
}


