#' Create a Venn Diagram from Long Data Format
#'
#' This function generates a Venn diagram based on a long format data frame containing at least
#' two columns: one representing the grouping variable and one representing the protein identifiers.
#' It uses the `eulerr` package to create Euler diagrams, which are similar to Venn diagrams, with
#' customizable colors and fill alpha values.
#'
#' @param dat_long A data frame in long format containing at least two columns:
#'   one representing the grouping variable (`group.by`), and one representing the
#'   protein identifiers (`proteincol.id`).
#' @param group.by A character string specifying the column name for the grouping variable
#'   (default is "sample"). This column will be used to divide the data into groups.
#' @param proteincol.id A character string specifying the column name for protein identifiers
#'   (default is "Accession"). This column will be used to identify unique proteins in the
#'   Venn diagram.
#' @param fill.alpha A numeric value between 0 and 1 that controls the transparency of the
#'   fill color in the Venn diagram (default is 0.3).
#' @param type plotting ellipse or circle(default is ellipse).
#'
#' @return A Venn diagram (Euler diagram) plot generated using the `eulerr` package. The plot
#'   shows the intersection of proteins across the different groups.
#'
#' @importFrom dplyr group_by distinct ungroup
#' @importFrom rlang sym
#' @import eulerr
#' @import ggplot2
#'
#' @examples
#' # Example
#' common_elements = c("a", "b", "c", "d", "e")
#' group1 = c(common_elements, sample(letters[6:26], 5))
#' group2 = c(common_elements, sample(letters[6:26], 5))
#' group3 = c(common_elements, sample(letters[6:26], 5))
#'
#' dt = data.frame(
#'   id = c(group1, group2, group3),
#'   group = rep(c("group1", "group2", "group3"), each = length(common_elements) + 5)
#' )
#' pltvenn(dt, group.by = "group", proteincol.id = "id",type=c("ellipse","circle")[1])
#'
#' @export
#'
pltvenn= function(dat_long, group.by = "sample", proteincol.id = "Accession",fill.alpha=0.3,type=c("ellipse","circle")[1]){

  #使用 sym() 和 !! 来动态引用列名
  group_col <- rlang::sym(group.by)
  protein_col <- rlang::sym(proteincol.id)

  #按指定的分组列去重
  proteins <- dat_long |>
    dplyr::group_by(!!group_col)  |>
    dplyr::distinct(!!protein_col, .keep_all = TRUE) |>
    dplyr::ungroup()

  #色
  my_fin_colors <-   c("turquoise3", "palevioletred3", "#ceca7c", "#c59fc9", "#84b59f", "cornflowerblue", "salmon3","#5698c4", "#d88c9a")

    # c("#BE507E", "#ED7D31", "#9183BA", "#C00000", "#29719D","#BE507E", "#ED7D31", "#9183BA", "#C00000", "#29719D")

  list_proteins <- split(proteins[[proteincol.id]], proteins[[group.by]])

  #绘
  proteins_venn <- plot(eulerr::euler(list_proteins, shape =type), # ellipse   circle
                        quantities = list(type = c("counts", "percent"),cex=0.8),
                        labels = list(font =2,col=my_fin_colors,cex=1.2),# font=1 2 3 4代表集合标签 普通文本 加粗文本 斜体文本 斜体加粗文本
                        fill = my_fin_colors,
                        edges=TRUE,#是否增加线条颜色
                        #edge_col= my_fin_colors ,
                        alpha =fill.alpha
                        );proteins_venn

  return(proteins_venn)
}
