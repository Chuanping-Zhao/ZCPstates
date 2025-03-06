#' Median Normalization of a Data Frame
#'
#' This function performs median normalization on a data frame by adjusting each numeric column so that 
#' its median becomes equal to the global median of all numeric column medians. The identifier column, 
#' specified by \code{id_col}, is excluded from the normalization.
#'
#' The normalization process involves the following steps:
#' \enumerate{
#'   \item Identify the numeric columns by excluding the specified identifier column.
#'   \item Compute the median of each numeric column, ignoring \code{NA} values.
#'   \item Calculate the global median of these column medians.
#'   \item Normalize each numeric column by dividing each value by its column's median and then multiplying by the global median.
#' }
#'
#' @param df A data frame containing the data to be normalized. It should include one identifier column and multiple numeric columns.
#' @param id_col A string specifying the name of the identifier column in \code{df} that should be excluded from normalization. Default is \code{"Glycopeptide"}.
#'
#' @return A data frame with the same structure as \code{df} where the numeric columns have been normalized so that 
#'   their medians equal the global median computed from the original data.
#'
#' @importFrom dplyr summarise
#' @importFrom dplyr across
#' @importFrom dplyr all_of
#' @importFrom dplyr mutate
#' @importFrom dplyr cur_column
#'
#' @examples
#' # demo
#' df <- data.frame(
#'   Glycopeptide = paste0("GP", 1:5),
#'   A = c(10, 20, 30, NA, 50),
#'   B = c(5, NA, 15, 25, 35),
#'   C = c(100, 200, 300, 400, 500)
#' )
#'
#' # Perform median normalization
#' df_norm <- median_normalization(df, id_col = "Glycopeptide")
#' print(df_norm)
#'
#' @export
#' 
median_normalization <- function(df, id_col = "Glycopeptide") {
  # 1. 确定数值列：除指定的 id 列以外的所有列
  numeric_cols <- setdiff(names(df), id_col)
  
  # 2. 计算各数值列的中位数（忽略 NA）
  col_medians <- df  |> 
    dplyr::summarise(dplyr::across(dplyr::all_of(numeric_cols), function(x) median(x, na.rm = TRUE))) |> 
    unlist()
  
  # 3. 计算所有列中位数的中位数，即全局中位数
  global_median <- median(col_medians, na.rm = TRUE)
  
  # 4. 对每个数值列进行中位数归一化
  df_norm <- df |> 
    dplyr::mutate( dplyr::across(
      dplyr::all_of(numeric_cols),
      ~ .x / col_medians[dplyr::cur_column()] * global_median
    ))
  
  return(df_norm)
}