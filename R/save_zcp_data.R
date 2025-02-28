
#' Batch Write Data Frames from a List to CSV Files
#'
#' This function takes a list of data frames and saves each data frame as a separate
#' CSV file in the specified directory. If the directory does not exist, it is created.
#' It uses the `purrr::walk` function to iterate over the list and save each `data.frame`
#' with its name as the CSV file name.
#'
#' @param dt_list A list where each element is a `data.frame`. The names of the list elements
#'   will be used as the names of the output CSV files.
#' @param save_dir A character string representing the directory path where the CSV files will be saved.
#'   If the directory does not exist, it will be created.
#'
#' @return This function does not return a value. It saves CSV files to the specified directory.
#'
#' @importFrom purrr walk
#'
#' @examples
#' # Example
#' df1 <- data.frame(a = 1:3, b = letters[1:3])
#' df2 <- data.frame(x = 4:6, y = letters[4:6])
#' batch_write_csv(dt_list=list(df1 = df1, df2 = df2), save_dir="outfile")
#'
#' @export
#'
batch_write_csv = function(dt_list, save_dir) {

  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }

  # 使用 walk 遍历 list，将每个 dataframe 保存到指定目录
  purrr::walk(names(dt_list), function(df_name) {
    df  = dt_list[[df_name]]
    # 定义输出文件路径
    output_file  = file.path(save_dir, paste0(df_name, ".csv"))

    # 保存 dataframe 到文件
    write.csv(df, output_file, row.names = FALSE)
    message("Saved: ", output_file)
  })
}
