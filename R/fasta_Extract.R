#' Extract Information from a FASTA File
#'
#' This function extracts relevant protein information such as UniProt ID, Gene Name, 
#' and Gene Description from a given FASTA file containing protein sequences.
#'
#' @param fasta.path A character string specifying the path to the FASTA file.
#' @return A tibble data frame with three columns:
#'   \item{Protein.Group}{UniProt ID extracted from the header.}
#'   \item{Genes}{Gene name extracted from the header.}
#'   \item{Description}{Gene description extracted from the header.}
#' @details
#' This function reads a FASTA file, extracts the UniProt ID, Gene Name, and Gene 
#' Description from each sequence header, and returns the extracted information in a tibble format.
#' The function uses regular expressions to extract these details based on the expected 
#' format of the FASTA headers.
#' 
#' The FASTA header is expected to follow the format:
#' \code{>sp|UniProt_ID|... GN=GeneName ... Description text OS=...}
#'
#' The function assumes that the header information contains at least the UniProt ID, 
#' Gene Name, and Gene Description for proper extraction.
#' 
#' @importFrom Biostrings readAAStringSet
#' @importFrom dplyr bind_rows
#' @importFrom purrr map
#' @importFrom stringr str_extract
#' @export
#'
#' @examples
#' # Example of extracting information from a FASTA file:
#' fasta_info <- Fasta_extract("your_fasta_file.fasta")
#' 
Fasta_extract=function(fasta.path="uniprotkb_proteome_UP000005640_AND_revi_2024_02_17.fasta"){
fasta=Biostrings::readAAStringSet(fasta.path)
#fasta.test=fasta[1]
fatsaLibs <- map(names(fasta), ~{
  names_info = .x
  #提取UniProt ID
  uniprot_id =  sub("^sp\\|([A-Za-z0-9]+)\\|.*", "\\1",names_info)
  #提取Gene Name
  gene_name = sub(".*GN=([A-Za-z0-9_]+).*", "\\1",names_info)
  #提取Gene Description   sub(".*_HUMAN\\s*(.*?)\\s*OS=.*", "\\1",names_info)
  gene_description =sub("^sp\\|[^|]*\\|[^ ]+\\s+(.*?)\\s+OS=.*", "\\1", names_info)  |> 
    trimws()
    
  tibble(Protein.Group = uniprot_id, Genes = gene_name, Description = gene_description)
  
})

#转换数据框
fatsaLibs <- bind_rows(fatsaLibs)

return(fatsaLibs)
}
