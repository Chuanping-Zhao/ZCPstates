#' @title uniprotEnrich: Functional Enrichment Analysis for Proteins
#' @description Perform pathway and GO enrichment analysis using custom protein datasets and libraries.
#' The function handles differential protein datasets, computes enrichment results, and saves output files.
#'
#' @param dt A data frame containing protein data. Must include the specified `protein.col.name` and `diff.condition.col.name` columns.
#' @param protein.col.name The column name in `dt` that contains protein IDs. Default is "ProteinIDs".
#' @param Protein.separator The separator used in the protein ID column. Set to `NULL` if no separator is used. Default is ";".
#' @param diff.condition.col.name The column name in `dt` that contains differential condition markers. Default is "sig".
#' @param diff.markers A character vector indicating the differential condition markers (e.g., `c("Up", "Down", "NotSig")`). Default is `c("Up", "Down", "NotSig")`.
#' @param library A data frame containing GO/Pathway annotations for proteins. Columns must include "Entry", "Pathway", and GO-related columns.
#' @param dt.include.allproteins Logical, whether `dt` includes all proteins (differential and non-differential). Default is `TRUE`.
#' @param cutoff Numeric, minimum number of proteins in a pathway/GO term to consider for enrichment. Default is `2`.
#' @param savefile The directory path where results will be saved. Default is "outputfile_test".
#' @return A list containing two elements for each condition (`Up`, `Down`, etc.):
#' \describe{
#'   \item{pathway_enrich_result}{A data frame with pathway enrichment results.}
#'   \item{go_enrichment}{A data frame with GO enrichment results (BP, CC, MF).}
#' }
#' @examples
#' # Example usage:
#' library(ZCPstates)
#' data("uniprotEnrich.demo.diff", package = "ZCPstates")
#' data("uniprotEnrich.demo.library.uniprot", package = "ZCPstates")
#' result <- uniprotEnrich(
#'   dt = uniprotEnrich.demo.diff,
#'   library=uniprotEnrich.demo.library.uniprot,
#'   protein.col.name="Entry",
#'   Protein.separator=";",
#'   diff.condition.col.name="Sig",
#'   diff.markers=c("Up","Down","NotSig"),
#'   dt.include.allproteins=TRUE,
#'    cutoff=1,
#'    savefile="outputfile"
#' )
#' @export
#' @importFrom purrr map_chr map walk2
#' @importFrom dplyr filter mutate select bind_rows rename everything
#' @importFrom stringr str_split str_trim
#' @importFrom tidyr separate_rows
#' @importFrom data.table as.data.table set setnames
#' @importFrom stats phyper p.adjust
#' @importFrom readr write_csv
#' @importFrom crayon yellow blue
#'
uniprotEnrich=function(
    dt,
    protein.col.name,
    Protein.separator=";",
    diff.condition.col.name,
    diff.markers=c("Up","Down","NotSig"),
    library,
    dt.include.allproteins=TRUE,
    cutoff=2,
    savefile="outputfile"

){

  if(!dir.exists(savefile)){
    dir.create(savefile)
  }else{
    cat(paste0("zcp: '",savefile,"' already exists!"))
  }



  required.cols= c(protein.col.name,diff.condition.col.name)
  missing.cols=required.cols[!required.cols %in% colnames(dt)]
  if (length(missing.cols) > 0) {
    stop( crayon::yellow(paste("The following columns are missing in 'dt':", paste(missing.cols, collapse = ", "),"--zcp")))
  }

  library.cols=colnames(library)
  lib.required.cols=c("Entry","Pathway","Gene Ontology (biological process)",
                      "Gene Ontology (cellular component)","Gene Ontology (molecular function)")
  lib.missing.cols = lib.required.cols[!lib.required.cols %in% library.cols]
  if (length(lib.missing.cols) > 0) {
    stop(crayon::yellow(paste("The library is missing the following columns:", paste(lib.missing.cols, collapse = ", "),"--zcp")))
  }

  if(cutoff>1000){warning( crayon::blue( paste0("The cutoff more than ",cutoff,",are you kidding me? --zcp")))}



  GO_dataset_dataframe=library#
  if(dt.include.allproteins==TRUE){

    if(!is.null(Protein.separator))  {
      background=unique(purrr::map_chr(stringr::str_split(dt[[protein.col.name]], Protein.separator), ~ .x[1]))

    }else{background=unique(dt[[protein.col.name]])}#

  }else{#
    background=GO_dataset_dataframe[["Entry"]]
  }
  #
  object_dataframe=GO_dataset_dataframe[which(GO_dataset_dataframe[["Entry"]]%in%background),]
  object_dataframe=object_dataframe[, c("Entry", "Pathway", "Gene Ontology (biological process)",
                                        "Gene Ontology (cellular component)", "Gene Ontology (molecular function)")]
  colnames(object_dataframe) <- c("Protein", "Pathway", "BP", "CC", "MF")


  #
  dt.processed = dt  |>
    dplyr::filter(get(diff.condition.col.name) != diff.markers[3])  |>
    dplyr::mutate(Protein = stringr::str_split(get(protein.col.name), ";") |>
                    purrr::map_chr(~ .x[1])) |>
    dplyr::filter(!is.na(Protein) & Protein != "NULL" & Protein != "Inf" & Protein != "0")


  dt.list=split(dt.processed,dt.processed[[diff.condition.col.name]] )
  dt.list$all=dt.processed#


  total_protein=length(background)


  .process_data <- function(dt.list, object_dataframe) {
    process_single <- function(subject_proteins) {
      #
      subject_all <- object_dataframe[which(object_dataframe$Protein %in% subject_proteins), ]
      subject_protein <- nrow(subject_all)#


      #
      .pathway_enrich <- function() {
        #

        object.database.clean=object_dataframe |>
          dplyr::select(Protein,Pathway)

        object.database.clean= data.table::as.data.table(object.database.clean)[
          !(is.na(Pathway) | Pathway == "" | is.null(Pathway) | is.finite(Pathway))
        ]#

        object.database.clean=object.database.clean |> tidyr::separate_rows(Pathway, sep = ";") #
        #
        data.table::set(object.database.clean, j = "Pathway", value = stringr::str_trim(object.database.clean$Pathway))#
        data.table::set(object.database.clean,j="Pathway",value = sub("^PATHWAY:\\s*", "",object.database.clean$Pathway))#
        data.table::set(object.database.clean,j = "Pathway", value = gsub("\\{[^}]*\\}", "", object.database.clean$Pathway)) #
        repeat {
          data.table::set(object.database.clean, j = "Pathway", value = gsub("\\.+\\s*$", "", object.database.clean$Pathway)) #
          #
          if (!any(grepl("\\.+\\s*$", object.database.clean$Pathway))) break
        }
        object.database.clean=unique(object.database.clean, by = c("Protein", "Pathway"))#


        #
        proteincounts.per.pahtway.all=data.table::as.data.table(object.database.clean)[, .N, by = Pathway]
        data.table::setnames(proteincounts.per.pahtway.all, "N", "proteincounts.per.pahtway.all")

        #
        diff.protein.database.clean=data.table::as.data.table(object.database.clean)[Protein %in% subject_proteins]
        proteincounts.per.pahtway.diff=data.table::as.data.table(diff.protein.database.clean)[, .N, by = Pathway]
        data.table::setnames(proteincounts.per.pahtway.diff, "N", "proteincounts.per.pahtway.diff")

        #
        proteincounts.merged <- data.table::as.data.table(proteincounts.per.pahtway.all)[
          data.table::as.data.table(proteincounts.per.pahtway.diff),
          on = "Pathway",
          nomatch = 0L,  #
          allow.cartesian = FALSE #
        ]
        #
        proteincounts.merged$background.proteincounts=total_protein
        #
        proteincounts.merged$diff.proteincounts=subject_protein

        result=data.table::as.data.table(proteincounts.merged)[, p_value := 1 - stats::phyper(
          proteincounts.per.pahtway.diff - 1,
          diff.proteincounts,
          background.proteincounts - diff.proteincounts,
          proteincounts.per.pahtway.all ,
          lower.tail=FALSE
        )]
        result=data.table::as.data.table(result)[, FDR := p.adjust(p_value, method = "BH")]
        result=data.table::as.data.table(result)[, Enrich_factor := round(proteincounts.per.pahtway.diff/proteincounts.per.pahtway.all,3)]
        data.table::setnames(result, old = "proteincounts.per.pahtway.diff", new = "counts")


        hits=data.table::as.data.table(diff.protein.database.clean)[, .(Protein = paste(Protein, collapse = "/")), by = Pathway]

        Final.results= data.table::as.data.table(result)[
          data.table::as.data.table(hits),
          on = "Pathway",
          nomatch = 0L,  #
          allow.cartesian = FALSE
        ]

        Final.results=data.table::as.data.table(Final.results)[, c("diff.proteincounts", "background.proteincounts") := NULL]
        data.table::setnames(Final.results, old = "Protein", new = "hited.proteins")
        Final.results$Type="Pathway Enrichment"
        return(as.data.frame(Final.results))
      }



      .GO_enrich=function(object_dataframe,Type){#Type="CC"

        object.database.clean=object_dataframe |>
          dplyr::select(Protein,any_of(Type)) |>
          dplyr::rename(Pathway = !!Type)

        object.database.clean= data.table::as.data.table(object.database.clean)[
          !(is.na(Pathway) | Pathway == "" | is.null(Pathway) | is.finite(Pathway))
        ]#

        object.database.clean=object.database.clean |> tidyr::separate_rows(Pathway, sep = ";") #
        #
        data.table::set(object.database.clean, j = "Pathway", value = stringr::str_trim(object.database.clean$Pathway))#

        repeat {
          data.table::set(object.database.clean, j = "Pathway", value = gsub("\\.+\\s*$", "", object.database.clean$Pathway)) #
          #
          if (!any(grepl("\\.+\\s*$", object.database.clean$Pathway))) break
        }
        object.database.clean=unique(object.database.clean, by = c("Protein", "Pathway"))#


        #
        proteincounts.per.pahtway.all=data.table::as.data.table(object.database.clean)[, .N, by = Pathway]
        data.table::setnames(proteincounts.per.pahtway.all, "N", "proteincounts.per.pahtway.all")

        #
        diff.protein.database.clean=data.table::as.data.table(object.database.clean)[Protein %in% subject_proteins]
        proteincounts.per.pahtway.diff=data.table::as.data.table(diff.protein.database.clean)[, .N, by = Pathway]
        data.table::setnames(proteincounts.per.pahtway.diff, "N", "proteincounts.per.pahtway.diff")

        #
        proteincounts.merged <- data.table::as.data.table(proteincounts.per.pahtway.all)[
          data.table::as.data.table(proteincounts.per.pahtway.diff),
          on = "Pathway",
          nomatch = 0L,  #
          allow.cartesian = FALSE #
        ]

        proteincounts.merged$background.proteincounts=total_protein

        proteincounts.merged$diff.proteincounts=subject_protein

        result=data.table::as.data.table(proteincounts.merged)[, p_value := 1 - stats::phyper(
          proteincounts.per.pahtway.diff - 1,
          diff.proteincounts,
          background.proteincounts - diff.proteincounts,
          proteincounts.per.pahtway.all ,
          lower.tail=FALSE
        )]
        result=data.table::as.data.table(result)[, FDR := p.adjust(p_value, method = "BH")]
        result=data.table::as.data.table(result)[, Enrich_factor := round(proteincounts.per.pahtway.diff/proteincounts.per.pahtway.all,3)]
        data.table::setnames(result, old = "proteincounts.per.pahtway.diff", new = "counts")


        hits=data.table::as.data.table(diff.protein.database.clean)[, .(Protein = paste(Protein, collapse = "/")), by = Pathway]#data.table::as.data.table(diff.protein.database.clean)[, .(Protein = paste(Protein, collapse = ";")), by = Pathway]

        Final.results= data.table::as.data.table(result)[
          data.table::as.data.table(hits),
          on = "Pathway",
          nomatch = 0L,
          allow.cartesian = FALSE #
        ]

        Final.results=data.table::as.data.table(Final.results)[, c("diff.proteincounts", "background.proteincounts") := NULL]
        data.table::setnames(Final.results, old = "Protein", new = "hited.proteins")
        Final.results$Type=Type

        return(Final.results)

      }

      pathway_enrich_result <- .pathway_enrich()
      bp_enriche_result <- .GO_enrich(object_dataframe =object_dataframe ,Type="BP")
      cc_enriche_result <- .GO_enrich(object_dataframe =object_dataframe ,Type="CC")
      mf_enriche_result <- .GO_enrich(object_dataframe =object_dataframe ,Type="MF")

      go_enrichment=dplyr::bind_rows(bp_enriche_result,cc_enriche_result,mf_enriche_result)



      return(list(
        pathway_enrich_result = pathway_enrich_result |>  dplyr::filter(counts>=cutoff),
        go_enrichment = go_enrichment|>  dplyr::filter(counts>=cutoff)
      ))
    }

    results <- purrr::map(dt.list, ~ process_single(.x$Protein))
    names(results) <- names(dt.list)
    return(results)
  }

  outputfile <- .process_data(dt.list, object_dataframe)



  purrr::walk2(
    .x = names(outputfile),
    .y = outputfile,
    .f = ~ {
      pathway_file <- file.path(savefile, paste0(.x, "_pathway_enrich_result.csv"))
      readr::write_csv(.y$pathway_enrich_result, pathway_file)
      go_file <- file.path(savefile, paste0(.x, "_go_enrich_result.csv"))
      readr::write_csv(.y$go_enrichment, go_file)
    }
  )

  cat( crayon::yellow("All results have been saved to", savefile, "--zcp\n"))

  return(outputfile)
}
