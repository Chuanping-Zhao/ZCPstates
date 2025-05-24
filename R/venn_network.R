#' Plot a Set-Based Venn Network with Side Tables
#'
#' This function visualizes overlapping sets as a network graph, with individual
#' elements and their memberships. Unique elements are also shown in side tables.
#'
#' @param df A data.frame or data.table.
#' @param protein_col Column name for protein IDs.
#' @param method_col Column name for method/group.
#' @param max_show Max proteins per side table.
#' @param font_size Font size for side tables.
#' @param font_face Font face (e.g., "plain", "bold").
#' @param font_family Font family for all text.
#' @param alpha Table background alpha (0–1).
#' @param show_labels Logical. Whether to show protein node labels.
#' @param show_table Logical. Whether to show unique element side tables.
#' @param set_colors Optional named vector of colors for each set (names must match methods).
#'
#' @return A patchwork plot with network and side tables (if `show_table = TRUE`), or network only.
#'
#' @examples
#' set.seed(123)
#' demo_proteins <- paste0("P", sprintf("%05d", sample(10000:99999, 600)))
#' group1 <- sample(demo_proteins, 300)
#' group2 <- sample(demo_proteins, 300)
#' df_demo <- data.frame(
#'   Protein = c(group1, group2),
#'   Method = c(rep("A", length(group1)), rep("B", length(group2)))
#' )
#'
#' venn_network(df_demo, show_labels = TRUE, show_table = TRUE)
#'
#' venn_network(df_demo, show_table = TRUE, set_colors = c(A = "#3A7A77", B = "#60C080"))
#'
#' @export
venn_network <- function(df, protein_col = "Protein", method_col = "Method", max_show = 10,
                         font_size = 10, font_face = "plain", font_family = "mono",
                         alpha = 0.3, show_labels = FALSE, show_table = FALSE,
                         set_colors = NULL) {
  
  make_colored_table <- function(df, bg_color, title = "Protein") {
    fill_colors <- rep(c("white", bg_color), length.out = nrow(df))
    rgba_fill <- scales::alpha(fill_colors, alpha)
    rgba_header <- scales::alpha(bg_color, min(alpha + 0.2, 1))
    
    colnames(df) <- title
    
    tg <- gridExtra::tableGrob(
      df, rows = NULL,
      theme = gridExtra::ttheme_default(
        core = list(
          bg_params = list(fill = rgba_fill, col = "grey60", lwd = 0.5),
          fg_params = list(
            col = "black", fontsize = font_size,
            fontface = font_face, fontfamily = font_family,
            hjust = 0, x = grid::unit(0.01, "npc")
          )
        ),
        colhead = list(
          bg_params = list(fill = rgba_header, col = NA),
          fg_params = list(
            col = "black", fontsize = font_size + 1,
            fontface = "bold", fontfamily = font_family
          )
        )
      )
    )
    patchwork::wrap_elements(full = tg)
  }
  
  data.table::setDT(df)
  df <- df[, .(Protein = get(protein_col), Method = get(method_col))]
  sets <- split(df$Protein, df$Method)
  set_names <- names(sets)
  
  if (is.null(set_colors)) {
    set_colors <- scales::hue_pal()(length(set_names))
    names(set_colors) <- set_names
  } else {
    if (!all(set_names %in% names(set_colors))) {
      stop("`set_colors` must be a named vector including all set names.")
    }
  }
  
  edges <- do.call(rbind, lapply(set_names, function(set) {
    data.frame(
      from = as.character(unname(sets[[set]])),
      to = rep(set, length(sets[[set]])),
      edge_color = rep(set_colors[set], length(sets[[set]])),
      stringsAsFactors = FALSE
    )
  }))
  
  elements <- unique(unlist(sets))
  protein_sets <- data.table::rbindlist(lapply(set_names, function(set) {
    data.table::data.table(Protein = sets[[set]], Set = set)
  }))
  protein_count <- protein_sets[, .(SetCount = .N, BelongsTo = paste(Set, collapse = ",")), by = Protein]
  
  nodes <- data.frame(
    name = c(set_names, elements),
    type = c(rep("Set", length(set_names)), rep("Element", length(elements))),
    stringsAsFactors = FALSE
  )
  nodes$size <- ifelse(nodes$type == "Set", 4, 1.5)
  nodes <- dplyr::left_join(nodes, protein_count, by = c("name" = "Protein"))
  
  nodes$group_layout <- dplyr::case_when(
    nodes$type == "Set" ~ nodes$name,
    nodes$type == "Element" & !is.na(nodes$SetCount) & nodes$SetCount == 1 ~ nodes$BelongsTo,
    TRUE ~ "Shared"
  )
  
  nodes$color <- "grey70"
  nodes$color[nodes$type == "Set"] <- set_colors[nodes$name[nodes$type == "Set"]]
  is_unique <- nodes$type == "Element" & !is.na(nodes$SetCount) & nodes$SetCount == 1
  nodes$color[is_unique] <- sapply(nodes$BelongsTo[is_unique], function(x) {
    set_colors[strsplit(x, ",")[[1]][1]]
  })
  

  g <- igraph::graph_from_data_frame(edges, vertices = nodes)
  layout <- ggraph::create_layout(g, layout = "kk")
  layout$group_layout <- nodes$group_layout[match(layout$name, nodes$name)]
  layout$y[layout$group_layout == set_names[1]] <- layout$y[layout$group_layout == set_names[1]] + 5
  layout$y[layout$group_layout == set_names[2]] <- layout$y[layout$group_layout == set_names[2]] - 5
  unique_df <- lapply(set_names, function(set) {
    unique(setdiff(sets[[set]], unlist(sets[set_names != set])))
  })
  unique_df <- lapply(unique_df, function(x) data.frame(Protein = head(x, max_show)))
  
  left_tbl <- make_colored_table(unique_df[[1]], set_colors[1], title = set_names[1])
  right_tbl <- make_colored_table(unique_df[[2]], set_colors[2], title = set_names[2])
  
  p_graph <- ggraph::ggraph(layout) +
    ggraph::geom_edge_link(ggplot2::aes(color = I(edge_color)), alpha = 0.1, width = 0.2) +
    ggraph::geom_node_point(ggplot2::aes(color = I(color), shape = type, size = size)) +
    ggraph::geom_node_text(data = dplyr::filter(layout, type == "Set"),
                           ggplot2::aes(label = name), size = 6, fontface = "bold", vjust = -1.5) +
    ggplot2::scale_size_identity() +
    ggplot2::scale_shape_manual(values = c(Set = 15, Element = 16)) +
    ggplot2::theme_void(base_family = font_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold", family = font_family),
      legend.position = "none",
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    ) +
    ggplot2::ggtitle("Protein Set Network")
  
  if (show_labels) {
    p_graph <- p_graph +
      ggraph::geom_node_text(
        data = dplyr::filter(layout, type == "Element"),
        ggplot2::aes(label = name),
        size = 2, hjust = -0.2, check_overlap = TRUE
      )
  }
  
  if (show_table) {
    left_tbl + p_graph + right_tbl + patchwork::plot_layout(widths = c(0.22, 0.56, 0.22))
  } else {
    p_graph
  }
}
