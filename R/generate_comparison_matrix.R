#' Generate All Pairwise Comparison Matrix
#'
#' This function creates a comparison matrix for all pairwise combinations of experimental conditions,
#' useful for differential analysis. It supports forcing specific groups to always appear in the denominator.
#'
#' @param sampleInfo A data frame or tibble containing at least one column specifying group information.
#' @param group_col Character. The name of the column in `sampleInfo` that contains the condition or group labels. Default is `"Condition"`.
#' @param force_denominator Optional character vector. Conditions specified here will always be placed in the denominator (i.e., assigned `-1` in the comparison vector) when involved in a comparison.
#'
#' @return A matrix with rows representing comparisons (e.g., `"A_vs_B"`) and columns corresponding to each unique condition. Each row is a contrast vector with `1` for numerator, `-1` for denominator, and `0` otherwise.
#'
#' @examples
#' sampleInfo <- tibble::tibble(
#'   Run = paste0("Run", 1:6),
#'   Condition = c("A", "A", "B", "B", "C", "C")
#' )
#' 
#' generate_comparison_matrix(sampleInfo)
#' generate_comparison_matrix(sampleInfo, force_denominator = c("C"))
#'
#' @export
#' 
generate_comparison_matrix <- function(sampleInfo, 
                                       group_col = "Condition", 
                                       force_denominator = NULL) {
  
  conditions = unique(sampleInfo[[group_col]])
  
  
  comb <- combn(conditions, 2, simplify = FALSE)
  
  comparison_list =list()
  
  for (pair in comb) {
    cond1 <- pair[1]
    cond2 <- pair[2]
    
    #默认 cond1 是分子+1，cond2是分母-1
    numerator = cond1
    denominator = cond2
    
    #如果其中某个是强制为分母的组，调整顺序
    if (!is.null(force_denominator)) {
      if (cond1 %in% force_denominator & !(cond2 %in% force_denominator)) {
        #使 cond1 作为分母
        numerator = cond2
        denominator =cond1
      } else if (cond2 %in% force_denominator & !(cond1 %in% force_denominator)) {
        #顺序正确，维持不变
      } else if (cond1 %in% force_denominator & cond2 %in% force_denominator) {
        #两者都在强制分母中，按默认顺序处理
      }
    }
    
    #创建向量
    vec = rep(0, length(conditions))
    names(vec) = conditions
    vec[numerator] =1
    vec[denominator]=-1
    
    comp_name = paste0(numerator, "_vs_", denominator)
    comparison_list[[comp_name]] =vec
  }
  
  #转换为矩阵
  comparison_matrix = do.call(rbind, comparison_list)
  row.names(comparison_matrix) = names(comparison_list)
  colnames(comparison_matrix) = conditions
  
  return(comparison_matrix)
}
