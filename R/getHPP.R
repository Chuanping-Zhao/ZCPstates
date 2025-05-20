#' @title Download and Parse Human Protein Atlas Blood Proteome Data
#'
#' @description
#' Downloads the Human Protein Atlas (HPA) blood proteome concentration dataset
#' (mass spectrometry or immunoassay based), unzips it, processes it into a tidy format,
#' computes log2 abundance, and assigns ranks based on abundance levels.
#'
#' @param dest_dir Character. Directory to save and extract the downloaded ZIP file.
#' Default is the current working directory (`getwd()`).
#'
#' @param technology Character. Choose between `"ms"` (mass spectrometry) or
#' `"immunoassay"` data. Default is `"ms"`.
#'
#' @param verbose Logical. Whether to print progress messages. Default is `TRUE`.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{`Data`}{A `data.table` containing the processed HPA data with columns:
#'     - `Gene`: Gene symbol
#'     - `ENSG ID`: Ensembl gene ID (only for `"ms"`)
#'     - `Conc [pg/L]`: Original concentration
#'     - `Technology`: `"ms"` or `"immunoassay"`
#'     - `Abundance`: log2-transformed abundance
#'     - `order`: Rank based on abundance
#'   }
#'   \item{`Data.version.Info`}{A character string with the data source URL.}
#' }
#'
#' @details
#' This function is useful for integrating HPA blood concentration data into proteomics workflows.
#' Internally uses `download.file()`, `unzip()`, and `data.table::fread()` for I/O operations.
#'
#' @source \url{https://www.proteinatlas.org/humanproteome/blood}
#'
#' @examples
#'
#' result = getHPP(dest_dir =getwd(), technology=c("immunoassay","ms")[2],verbose = TRUE)
#' head(result$Data)
#'
#' @export
#'
getHPP = function(dest_dir =getwd(), technology=c("immunoassay","ms")[2],verbose = TRUE) {
  technology =match.arg(technology, choices = c("ms", "immunoassay"))
  #url
  url =switch(technology,
                ms = "https://www.proteinatlas.org/download/tsv/blood_ms_concentration.tsv.zip",
                immunoassay = "https://www.proteinatlas.org/download/tsv/blood_immunoassay_concentration.tsv.zip")

  zip_name = paste0("blood_", technology, "_concentration.tsv.zip")
  zip_path = file.path(dest_dir, zip_name)


  if (verbose) message("Downloading data from HPA...--zcp")
  tryCatch({
    invisible(capture.output(
      download.file(url, destfile = zip_path, mode = "wb", quiet = TRUE)
    ))
  }, error = function(e) {
    stop("Download failed. Please check your internet connection or the URL.--zcp")
  })

  if (verbose) message("📦 Unzipping file...--zcp")
  unzip(zip_path, exdir = dest_dir)
  tsv_file = list.files(dest_dir, pattern = "\\.tsv$", full.names = TRUE)
  if (length(tsv_file) == 0) stop("No .tsv file found after unzip.--zcp")

  if(technology=="ms"){
  HPP =data.table::fread(tsv_file[1])
  HPP[,Technology:=technology]
  conc_col= grep("^Conc", names(HPP), value = TRUE,ignore.case = TRUE)
  HPP[, Abundance := log2(as.numeric(get(conc_col)))]
  data.table::setnames(HPP, old = conc_col, new = "Conc_pg_per_L")
  HPP = HPP[order(-Abundance)][, order := .I]#  "Gene","ENSG ID","Conc [pg/L]", "Technology", "Abundance","order"
  }else{
    HPP =data.table::fread(tsv_file[1])
    data.table::setnames(x = HPP,old = c("gene"),new = c("Gene"))
    HPP[,Technology:=technology]
    conc_col= grep("^Conc", names(HPP), value = TRUE,ignore.case = TRUE)
    HPP[, Abundance := log2(as.numeric(get(conc_col)))]
    data.table::setnames(HPP, old = conc_col, new = "Conc_pg_per_L")
    HPP = HPP[order(-Abundance)][, order := .I]
  }

  suppressWarnings({
    file.remove(zip_path)
    file.remove(tsv_file)
  })
  #load("HumanFastaLibs.rda")
  data("HumanFastaLibs", envir = environment())
  data.table::setDT(HPP)
  data.table::setDT(HumanFastaLibs)
  HPP_annotated= merge(HPP, HumanFastaLibs, by = "Gene", all.x = TRUE)


  dataInfo="Source: https://www.proteinatlas.org/humanproteome/blood"
  if (verbose) {
    message("✅ Done. ", nrow(HPP), " rows and ", ncol(HPP), " columns loaded.")
  }

  return(list(Data=HPP_annotated,Data.version.Info=dataInfo))
}

