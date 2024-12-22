#' @title Save Plots in Multiple Formats
#' @description This function saves a given ggplot object in PDF, PNG, and optionally PPT formats.
#' It ensures the output directory exists and generates the files in the specified directory.
#' @param Fig A `ggplot` object to save.
#' @param FigName A string representing the base name for the saved files (without extension).
#' @param outputfile A string representing the directory where the files will be saved. The function creates the directory if it does not exist.
#' @param widths A numeric value specifying the width of the saved plot.
#' @param heights A numeric value specifying the height of the saved plot.
#' @param ppt Logical. If `TRUE` (default), the function generates a PPT file in addition to PNG and PDF files. If `FALSE`, only PNG and PDF files are created.
#' @return No return value. The function saves the plot files in the specified directory.
#' @examples
#' \dontrun{
#' library(ggplot2)
#' #demo
#' example_plot <- ggplot(mtcars, aes(mpg, wt)) + geom_point()
#'
#' # Save the plot
#' save_zcp(Fig = example_plot,
#'          FigName = "example_plot",
#'          outputfile = "outputfile",
#'          widths = 7,
#'          heights = 5,
#'          ppt = TRUE)
#' }
#' @importFrom ggplot2 ggsave
#' @importFrom export graph2ppt
#' @export
save_zcp <- function(Fig, FigName, outputfile, widths, heights, ppt = TRUE) {
  if (dir.exists(outputfile)) {

    cat(paste("zcp:", outputfile, "already exists!\n"))
  } else {

    dir.create(outputfile, recursive = TRUE)
    cat(paste("zcp:", outputfile, "has been created!\n"))
  }


  Filepaths <- paste0(outputfile, "/", FigName, c(".pdf", ".png"))
  ggplot2::ggsave(Filepaths[1], width = widths, plot = Fig, height = heights, device = "pdf")
  ggplot2::ggsave(Filepaths[2], width = widths, plot = Fig, height = heights, device = "png")


  if (ppt) {
    ppt_filepath <- paste0(outputfile, "/", FigName, ".ppt")
    export::graph2ppt(x = Fig, file = ppt_filepath,
                      vector.graphic = TRUE, width = widths, height = heights, aspectr = sqrt(2), append = TRUE)
    cat(paste("zcp: PPT file generated at", ppt_filepath, "\n"))
  }
}
