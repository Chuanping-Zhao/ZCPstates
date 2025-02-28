#' Create a Directory if it Doesn't Exist
#'
#' This function checks if the specified directory exists. If it does not exist, 
#' it creates the directory. If no directory name is provided, it defaults to 
#' creating a directory named "CleanedData".
#'
#' @param x A character string representing the directory path. If `NULL`, the 
#'   function defaults to creating a directory named "CleanedData".
#' 
#' @return NULL
#' 
#' @examples
#' # Create a directory "CleanedData" (default)
#' dir_c()
#' 
#' # Create a directory "testFile"
#' dir_c("testFile")
#' @export
#'
dir_c=function(x = NULL){

  if (is.null(x)) {
    x= "CleanedData"
  }
  
  if (dir.exists(x)) {
    
    cat(paste("zcp:", x, "already exists!\n"))
  } else {
    
    dir.create(x, recursive = TRUE)
    cat(paste("zcp:", x, "has been created!\n"))
  }
  
  }
