#' easyOrgdbEnrich: An Easy Wrapper for Functional Enrichment & GSEA Analysis
#'
#' This function performs functional enrichment analysis (KEGG & GO) and GSEA analysis
#' on proteomics data using UniProt IDs. It supports organism-specific databases (`org.*.eg.db`)
#' and provides automatic ID conversion, annotation, visualization (bubble plots and GSEA plots),
#' and result export.
#'
#' @param dt A data.frame or tibble. Input differential expression results. Must include columns: `Protein`, `Label`, `log2FC`.
#' @param proteincol.id Character. Column name in `dt` containing UniProt protein IDs. Default is `"Protein"`.
#' @param comparison Character. Column name in `dt` representing comparison groups. Default is `"Label"`.
#' @param savepath Character. Path to save enrichment results and plots. Default is `"Report/enrich"`.
#' @param sep.pattern Character. Regex pattern to split multi-ID strings. Default is `";"`.
#' @param species Character. OrgDb database name (e.g., `"org.Hs.eg.db"`). Default is `"org.Hs.eg.db"`.
#' @param background A background dataframe with columns: `Protein`, `Label`, and `log2FC`, used as universe for enrichment. Default is NULL.
#' @param log2fc.col Character. Column name representing log2 fold-change values. Default is `"log2FC"`.
#' @param gsea_widths Numeric. Width of GSEA plots. Default is `4`.
#' @param gsea_heights Numeric. Height of GSEA plots. Default is `3`.
#' @param plotppt Logical. Whether to save figures as PPT using `save_zcp()`. Default is `FALSE`.
#'
#' @return No value is returned. The function saves enrichment results and GSEA figures to `savepath`.
#'
#' @details
#' The function includes:
#' \itemize{
#'   \item KEGG enrichment analysis using `clusterProfiler::enrichKEGG()`
#'   \item GO enrichment analysis using `clusterProfiler::enrichGO()`
#'   \item GSEA analysis using `clusterProfiler::gseGO()` and `gseKEGG()`
#'   \item Protein ID conversion via `clusterProfiler::bitr()`
#'   \item Bubble plots using `ggplot2`
#'   \item GSEA curve plotting using `enrichplot::gseaplot2()`
#'   \item Automatic saving using `ZCPstates::save_zcp()`
#' }
#'
#' @importFrom dplyr filter distinct select pull inner_join arrange group_by mutate desc
#' @importFrom tidyr separate
#' @importFrom tibble tibble deframe
#' @importFrom clusterProfiler enrichGO enrichKEGG gseGO gseKEGG bitr
#' @importFrom enrichplot gseaplot2
#' @importFrom ggplot2 ggplot aes geom_point theme_bw xlab ylab scale_color_gradient
#'             ggtitle theme element_text element_blank scale_size guide_colorbar
#' @importFrom purrr walk2 imap
#' @importFrom readr write_csv
#'
#' @examples
#' \dontrun{
#' easyOrgdbEnrich(
#'   dt = diff_sig,
#'   proteincol.id = "Protein",
#'   comparison = "Label",
#'   savepath = "Report/enrich",
#'   sep.pattern = ";",
#'   species = "org.Hs.eg.db",
#'   background = NULL,
#'   log2fc.col = "log2FC"
#' )
#' }
#'
#' @export

easyOrgdbEnrich = function(
    dt,# dt=diff_sig
    proteincol.id = "Protein",
    comparison = "Label",
    savepath = "Report/enrich",
    sep.pattern = ";",
    species = c("org.Hs.eg.db","org.Mm.eg.db","org.Rn.eg.db","org.Dm.eg.db","org.Ce.eg.db","org.Sc.sgd.db",
                "org.Dr.eg.db", "org.Bt.eg.db","org.Gg.eg.db","org.Cf.eg.db","org.Ss.eg.db")[1],
    background = NULL  ,
    log2fc.col = "log2FC",
    gsea_widths = 4,
    gsea_heights = 3,
    plotppt = FALSE
) {


  .check_species <- function(species) {
    valid_species <- c(
      "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "org.Dm.eg.db", "org.Ce.eg.db",
      "org.Sc.sgd.db", "org.Dr.eg.db", "org.Bt.eg.db", "org.Gg.eg.db", "org.Cf.eg.db", "org.Ss.eg.db"
    )
    if (!(species %in% valid_species)) {
      stop(paste0("Unsupported species: ", species))
    }
  }
  .check_species(species)


  dt.list = split(dt[[proteincol.id]], dt[[comparison]]) |>
    lapply(function(x) data.frame(proteins = unique(x))) |>
    (\(x) Filter(function(x) nrow(x) > 1, x))()


  id.trans = lapply(dt.list, function(x) {
    x |>
      dplyr::select(protein_groups = proteins) |>
      tidyr::separate(protein_groups, into = c("proteins", 'temp'), sep = sep.pattern) |>
      dplyr::select(proteins) |>
      dplyr::pull() |>
      clusterProfiler::bitr(fromType = "UNIPROT",
                            toType = c("ENSEMBL", "ENTREZID"),
                            OrgDb = species) |>
      dplyr::distinct(UNIPROT, .keep_all = TRUE)
  })


  names(id.trans) = names(dt.list)
  dirs = file.path(savepath, names(id.trans))
  dirs_figures_gseaGO= file.path(savepath, names(id.trans),"Figure_gseaGO")
  dirs_figures_gseaKEGG= file.path(savepath, names(id.trans),"Figure_gseaKEGG")
  .ensure_dirs <- function(dirs, recursive = TRUE, verbose = TRUE) {
    if (!is.character(dirs)) stop("'dirs' must be a character vector--zcp")
    new_dirs <- dirs[!dir.exists(dirs)]
    if (length(new_dirs) > 0) {
      for (dir in new_dirs) {
        success <- tryCatch({
          dir.create(dir, recursive = recursive, showWarnings = FALSE)
          TRUE
        }, error = function(e) {
          message("Error creating directory: ", dir, " Error: ", e$message)
          FALSE
        })
        if (success && verbose) {
          message("Created directory: ", dir)
        }
      }
    } else {
      if (verbose) message("All directories already exist.--zcp")
    }
    invisible(TRUE)
  }

  .ensure_dirs(dirs)
  .ensure_dirs(dirs_figures_gseaGO)
  .ensure_dirs(dirs_figures_gseaKEGG)


  .load_species_db <- function(species) {
    if (!require(species, character.only = TRUE)) {
      stop("Failed to load species package: ", species)
    }
    get(species)
  }


  kegg.filenames.id = paste0(savepath, "/", names(id.trans), "/kegg_id_mapping.csv")
  purrr::walk2(id.trans, kegg.filenames.id, readr::write_csv)


  .get_kegg_organism = function(species) {
    mapping <- c(
      "org.Hs.eg.db" = "hsa", "org.Mm.eg.db" = "mmu", "org.Rn.eg.db" = "rno", "org.Dm.eg.db" = "dme",
      "org.Ce.eg.db" = "cel", "org.Sc.sgd.db" = "sce", "org.Dr.eg.db" = "dre",
      "org.Bt.eg.db" = "bta", "org.Gg.eg.db" = "gga", "org.Cf.eg.db" = "cfa", "org.Ss.eg.db" = "ssc"
    )
    if (!(species %in% names(mapping))) stop("Unsupported species: ", species)
    return(mapping[species])
  }
  species_kegg = .get_kegg_organism(species)


  background.list <- list()
  if (!is.null(background)) {
    # Check required columns
    required_cols <- c("Protein", "Label", "log2FC")
    if (!all(required_cols %in% colnames(background))) {
      stop("If `background` is provided, it must contain columns: Protein, Label, log2FC.")
    }


    bg.split = split(background$Protein, background$Label)
    background.list = lapply(bg.split, function(x) {
      tibble::tibble(protein_groups = unique(x)) |>
        tidyr::separate(protein_groups, into = c("proteins", "temp"), sep = sep.pattern) |>
        dplyr::select(proteins) |>
        dplyr::pull() |>
        clusterProfiler::bitr(fromType = "UNIPROT",
                              toType = c("ENSEMBL", "ENTREZID"),
                              OrgDb = species) |>
        dplyr::distinct(UNIPROT, .keep_all = TRUE)
    })
  }


  kegg.list = list()
  if (is.null(getOption("clusterProfiler.download.method"))) {
    options(clusterProfiler.download.method = "auto")
  }
  options(timeout = max(600, getOption("timeout")))


  for (y in seq_along(id.trans)) {
    this_name = names(id.trans)[y]
    gene.df = id.trans[[y]]
    universe.df = if (!is.null(background)) background.list[[this_name]] else NULL
    result_df = tryCatch({
      clusterProfiler::enrichKEGG(
        gene = gene.df$ENTREZID,
        universe = if (!is.null(universe.df)) universe.df$ENTREZID else NULL,
        keyType = "kegg",
        organism = species_kegg,
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      )@result
    }, error = function(e) NULL)

    if (!is.null(result_df) && nrow(result_df) > 0) {
      KEGGData =result_df |> dplyr::arrange(dplyr::desc(Count)) |> dplyr::slice_head(n = 10)
      kegg_plot_top10 = ggplot2::ggplot(KEGGData,  ggplot2::aes(Count, Description, size = Count, color = pvalue)) +
        ggplot2::geom_point() +
        ggplot2::theme_bw() +
        ggplot2::xlab("Counts") +
        ggplot2::ylab(NULL) +
        ggplot2::scale_color_gradient(low = "#B61B47", high = "#4272B4",guide = ggplot2::guide_colorbar(
          barwidth = 1,   #宽度
          barheight = 4     #高度
        ))+
        ggplot2:: ggtitle(paste("KEGG analysis -", this_name)) +
        ggplot2::theme(plot.title =  ggplot2::element_text(size = 10, hjust = 0.5),
              plot.subtitle =  ggplot2::element_text(size = 10, hjust = 0.5),
              panel.grid.major.y=  ggplot2::element_blank(),
              panel.grid.minor.y= ggplot2::element_blank(),  #隐藏y轴次要网格线
              panel.grid.major.x= ggplot2::element_blank(),  #隐藏x轴主要网格线
              panel.grid.minor.x= ggplot2::element_blank()) +
        ggplot2::scale_size(range = c(3, 5))

      kegg.filenames.plot =base::file.path(savepath, this_name)
     save_zcp(Fig = kegg_plot_top10,FigName = "KEGG_top10_by_counts",outputfile =kegg.filenames.plot,widths = 4*1.5,heights = 4,ppt =plotppt  )
    }
   # result_df = enrich_obj@result


    id_map <- gene.df[, c("ENTREZID", "UNIPROT")] |> dplyr::distinct()
    result_df$geneID <- sapply(strsplit(result_df$geneID, "/"), function(ids) {
      mapped <- id_map$UNIPROT[match(ids, id_map$ENTREZID)]
      mapped <- mapped[!is.na(mapped)]
      paste(mapped, collapse = "/")
    })

    kegg.list[[y]] <- result_df
  }
  names(kegg.list) = names(id.trans)


  kegg.filenames = paste0(savepath, "/", names(id.trans), "/KEGG.csv")
  purrr::walk2(kegg.list, kegg.filenames, readr::write_csv)


  #go

    go.list = list()

  for (y in seq_along(id.trans)) {
    this_name = names(id.trans)[y]
    gene.df = id.trans[[y]]
    universe.df = if (!is.null(background)) background.list[[this_name]] else NULL

    result_df = tryCatch({
      clusterProfiler::enrichGO(
        gene = gene.df[["ENTREZID"]],
        #universe = if (!is.null(universe.df)) universe.df$ENTREZID else NULL,
        OrgDb = species,#.load_species_db(species)
        keyType = "ENTREZID",
        ont = "ALL",
        pAdjustMethod = 'fdr',
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = FALSE
      )@result
    }, error = function(e) NULL)

    #可视化
    if (!is.null(result_df) && nrow(result_df) > 0) {
      GOData =result_df |>dplyr::group_by(ONTOLOGY)  |> dplyr::arrange(dplyr::desc(Count)) |> dplyr::slice_head(n = 10)
      GO_plot_top10_BP = ggplot2::ggplot(GOData |> dplyr::filter(ONTOLOGY=="BP"),  ggplot2::aes(Count, Description, size = Count, color = pvalue)) +
        ggplot2::geom_point() +
        ggplot2::theme_bw() +
        ggplot2::xlab("Counts") +
        ggplot2::ylab(NULL) +
        ggplot2::scale_color_gradient(low = "#B61B47", high = "#4272B4",guide = ggplot2::guide_colorbar(
          barwidth = 1,   #宽度
          barheight = 4     #高度
        )) +
        ggplot2:: ggtitle(paste("GO analysis[BP] -", this_name)) +
        ggplot2::theme(plot.title =  ggplot2::element_text(size = 10, hjust = 0.5),
                       plot.subtitle =  ggplot2::element_text(size = 10, hjust = 0.5),
                       panel.grid.major.y=  ggplot2::element_blank(),
                       panel.grid.minor.y= ggplot2::element_blank(),  #隐藏y轴次要网格线
                       panel.grid.major.x= ggplot2::element_blank(),  #隐藏x轴主要网格线
                       panel.grid.minor.x= ggplot2::element_blank()) +
        ggplot2::scale_size(range = c(3, 5))

      GO_plot_top10_CC = ggplot2::ggplot(GOData |> dplyr::filter(ONTOLOGY=="CC"),  ggplot2::aes(Count, Description, size = Count, color = pvalue)) +
        ggplot2::geom_point() +
        ggplot2::theme_bw() +
        ggplot2::xlab("Counts") +
        ggplot2::ylab(NULL) +
        ggplot2::scale_color_gradient(low = "#B61B47", high = "#4272B4",guide = ggplot2::guide_colorbar(
          barwidth = 1,   #宽度
          barheight = 4     #高度
        )) +
        ggplot2:: ggtitle(paste("GO analysis[CC] -", this_name)) +
        ggplot2::theme(plot.title =  ggplot2::element_text(size = 10, hjust = 0.5),
                       plot.subtitle =  ggplot2::element_text(size = 10, hjust = 0.5),
                       panel.grid.major.y=  ggplot2::element_blank(),
                       panel.grid.minor.y= ggplot2::element_blank(),  #隐藏y轴次要网格线
                       panel.grid.major.x= ggplot2::element_blank(),  #隐藏x轴主要网格线
                       panel.grid.minor.x= ggplot2::element_blank()) +
        ggplot2::scale_size(range = c(3, 5))

      GO_plot_top10_MF = ggplot2::ggplot(GOData |> dplyr::filter(ONTOLOGY=="MF"),  ggplot2::aes(Count, Description, size = Count, color = pvalue)) +
        ggplot2::geom_point() +
        ggplot2::theme_bw() +
        ggplot2::xlab("Counts") +
        ggplot2::ylab(NULL) +
        ggplot2::scale_color_gradient(low = "#B61B47", high = "#4272B4",guide = ggplot2::guide_colorbar(
          barwidth = 1,   #宽度
          barheight = 4     #高度
        )) +
        ggplot2:: ggtitle(paste("GO analysis[MF] -", this_name)) +
        ggplot2::theme(plot.title =  ggplot2::element_text(size = 10, hjust = 0.5),
                       plot.subtitle =  ggplot2::element_text(size = 10, hjust = 0.5),
                       panel.grid.major.y=  ggplot2::element_blank(),
                       panel.grid.minor.y= ggplot2::element_blank(),  #隐藏y轴次要网格线
                       panel.grid.major.x= ggplot2::element_blank(),  #隐藏x轴主要网格线
                       panel.grid.minor.x= ggplot2::element_blank()) +
        ggplot2::scale_size(range = c(3, 5))

      kegg.filenames.plot =base::file.path(savepath, this_name)
      save_zcp(Fig = GO_plot_top10_BP,FigName = "GOBP_top10_by_counts",outputfile =kegg.filenames.plot,widths = 4*1.5,heights = 4,ppt =plotppt  )
      save_zcp(Fig = GO_plot_top10_CC,FigName = "GOCC_top10_by_counts",outputfile =kegg.filenames.plot,widths = 4*1.5,heights = 4,ppt =plotppt  )
      save_zcp(Fig = GO_plot_top10_MF,FigName = "GOMF_top10_by_counts",outputfile =kegg.filenames.plot,widths = 4*1.5,heights = 4,ppt =plotppt  )
    }




    # --- Convert geneID from ENTREZID to UNIPROT ---
    id_map <- gene.df[, c("ENTREZID", "UNIPROT")] |> dplyr::distinct()
    result_df$geneID <- sapply(strsplit(result_df$geneID, "/"), function(ids) {
      mapped <- id_map$UNIPROT[match(ids, id_map$ENTREZID)]
      mapped <- mapped[!is.na(mapped)]
      paste(mapped, collapse = "/")
    })

    go.list[[y]] <- result_df
  }
  names(go.list) = names(id.trans)

  # --- Save GO results ---
  go.filenames = paste0(savepath, "/", names(id.trans), "/GO.csv")
  purrr::walk2(go.list, go.filenames, readr::write_csv)


  #gsea_go
  gsea.list = list()
  gsea.visual.go=list()
  for (y in seq_along(id.trans)) {
    this_name = names(id.trans)[y]

logfc.vector=dt |>
  dplyr::filter(.data[[comparison]] == this_name) |>
  tidyr::separate(.data[[proteincol.id]], into = c("proteins", "temp"), sep = sep.pattern) |>
  dplyr::distinct(proteins, .keep_all = TRUE)

    geneList = clusterProfiler::bitr(logfc.vector$proteins,
                                     fromType = "UNIPROT",
                                     toType =c("ENSEMBL", "ENTREZID"),
                                     OrgDb = species) |>
      dplyr::inner_join(logfc.vector, by = c("UNIPROT" = "proteins")) |>
      dplyr::select(ENTREZID, !!rlang::sym(log2fc.col)) |>
      dplyr::distinct(ENTREZID, .keep_all = TRUE) |>
      dplyr::arrange(desc(!!rlang::sym(log2fc.col))) |>
      tibble::deframe()

    if (length(geneList) == 0) {
      message("GSEA: No geneList for group ", this_name)
      gsea.list[[y]] = NULL
      next
    }
    gsea.kegg.go=tryCatch({
      clusterProfiler::gseGO(
        geneList = geneList,
        OrgDb = .load_species_db(species),
        keyType = "ENTREZID",
        ont = "ALL",
        pvalueCutoff = 0.05,
        verbose = FALSE,
        minGSSize = 2
      )
    }, error = function(e) {
      message("GSEA error for group ", this_name, ": ", e$message)
      NULL
    })

    gsea.res = tryCatch({
      gsea.kegg.go@result
    }, error = function(e) {
      message("GSEA error for group ", this_name, ": ", e$message)
      NULL
    })


    if (!is.null(gsea.res) && nrow(gsea.res) > 0 && "core_enrichment" %in% colnames(gsea.res)) {
      id_map = clusterProfiler::bitr(logfc.vector$proteins,
                                     fromType = "UNIPROT",
                                     toType = "ENTREZID",
                                     OrgDb = species) |>
        dplyr::distinct(ENTREZID, .keep_all = TRUE)
      gsea.res$core_enrichment <- sapply(strsplit(gsea.res$core_enrichment, "/"), function(ids) {
        mapped <- id_map$UNIPROT[match(ids, id_map$ENTREZID)]
        paste(na.omit(mapped), collapse = "/")
      })
    }

    gsea.list[[y]] = gsea.res
    gsea.visual.go[[y]]=gsea.kegg.go
  }
  names(gsea.list) = names(id.trans)
  names(gsea.visual.go) = names(id.trans)

  gsea.filenames = paste0(savepath, "/", names(id.trans), "/GSEA_GO.csv")
  purrr::walk2(gsea.list, gsea.filenames, function(x, f) {
    if (!is.null(x)) readr::write_csv(x, f)
  })
  #批量可视化
  purrr::imap(gsea.visual.go, function(obj, name) {
    if (!is.null(obj) && nrow(obj@result) > 0) {

      output_dir <- base::file.path(savepath, name,"Figure_gseaGO")

      for (i in seq_len(nrow(obj@result))) {
        pathway_name <- obj@result$Description[i]
        clean_name <- gsub("[^A-Za-z0-9_]+", "_", pathway_name)  #合法化文件名
        fig_name <- paste0("GSEA_GO_", sprintf("%02d", i), "_", clean_name)

        p <- enrichplot::gseaplot2(obj, geneSetID = i, title = pathway_name)

        save_zcp(
          Fig = p,
          FigName = fig_name,
          outputfile = output_dir,
          widths = gsea_widths,
          heights = gsea_heights,
          ppt = plotppt
        )

        base::message("Saved plot for ", name, " → ", fig_name)
      }

    } else {
      base::message("Skipped ", name, " — empty or NULL gseaResult --zcp")
    }
  })




  #gsea_kegg


  gsea.list.kegg = list()
  gsea.visual.kegg=list()
  for (y in seq_along(id.trans)) {
    this_name = names(id.trans)[y]



    logfc.vector=dt |>
      dplyr::filter(.data[[comparison]] == this_name) |>
      tidyr::separate(.data[[proteincol.id]], into = c("proteins", "temp"), sep = sep.pattern) |>
      dplyr::distinct(proteins, .keep_all = TRUE)

    geneList = clusterProfiler::bitr(logfc.vector$proteins,
                                     fromType = "UNIPROT",
                                     toType =c("ENSEMBL", "ENTREZID"),
                                     OrgDb = species) |>
      dplyr::inner_join(logfc.vector, by = c("UNIPROT" = "proteins")) |>
      dplyr::select(ENTREZID, !!rlang::sym(log2fc.col)) |>
      dplyr::distinct(ENTREZID, .keep_all = TRUE) |>
      dplyr::arrange(desc(!!rlang::sym(log2fc.col))) |>
      tibble::deframe()

    if (length(geneList) == 0) {
      message("GSEA: No geneList for group ", this_name)
      gsea.list[[y]] = NULL
      next
    }

    gsea.kegg=tryCatch({
      clusterProfiler::gseKEGG(
        geneList = geneList,
        organism = species_kegg,
        keyType = "kegg",
        pvalueCutoff = 0.05,
        nPerm = 1000, # 置换次数
        minGSSize = 2, #基因集最小尺寸
        verbose = FALSE
      )
    }, error = function(e) {
      message("gseKEGG error for group ", this_name, ": ", e$message)
      NULL
    })

    gsea.res = tryCatch({gsea.kegg@result
    }, error = function(e) {
      message("gseKEGG error for group ", this_name, ": ", e$message)
      NULL
    })



    if (!is.null(gsea.res) && nrow(gsea.res) > 0 && "core_enrichment" %in% colnames(gsea.res)) {
      id_map = clusterProfiler::bitr(logfc.vector$proteins,
                                     fromType = "UNIPROT",
                                     toType = "ENTREZID",
                                     OrgDb = species) |>
        dplyr::distinct(ENTREZID, .keep_all = TRUE)
      gsea.res$core_enrichment <- sapply(strsplit(gsea.res$core_enrichment, "/"), function(ids) {
        mapped <- id_map$UNIPROT[match(ids, id_map$ENTREZID)]
        paste(na.omit(mapped), collapse = "/")
      })
    }

    gsea.list.kegg[[y]] = gsea.res
    gsea.visual.kegg[[y]]=gsea.kegg
  }
  names(gsea.list.kegg) = names(id.trans)
  names(gsea.visual.kegg) = names(id.trans)

  gsea.filenames = paste0(savepath, "/", names(id.trans), "/GSEA_KEGG.csv")
  purrr::walk2(gsea.list.kegg, gsea.filenames, function(x, f) {
    if (!is.null(x)) readr::write_csv(x, f)
  })


  #批量绘制gsea_kegg图形
  purrr::imap(gsea.visual.kegg, function(obj, name) {
    if (!is.null(obj) && nrow(obj@result) > 0) {

      output_dir <- base::file.path(savepath, name,"Figure_gseaKEGG")

      for (i in seq_len(nrow(obj@result))) {
        pathway_name <- obj@result$Description[i]
        clean_name <- gsub("[^A-Za-z0-9_]+", "_", pathway_name)  #合法化文件名
        fig_name <- paste0("GSEA_GO_", sprintf("%02d", i), "_", clean_name)

        p <- enrichplot::gseaplot2(obj, geneSetID = i, title = pathway_name)

       save_zcp(
          Fig = p,
          FigName = fig_name,
          outputfile = output_dir,
          widths = gsea_widths,
          heights = gsea_heights,
          ppt = plotppt
        )

        base::message("Saved plot for ", name, " → ", fig_name)
      }

    } else {
      base::message("Skipped ", name, " — empty or NULL gseaResult --zcp")
    }
  })



}



