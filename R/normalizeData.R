#' Data Normalization
#'
#' Normalize data to have mean zero and std 1
#'
#' @param dataset The input list of data sets matrix
#' @param enable_normalization An argument to decide whether to use normalizaiton or not,  default is TRUE
#' @param column_sum_normalization An argument to decide whether to use column sum normalization or not, default it TRUE
#' @param nonnegative_normalization An argument to decide wehther it is nonnegative matrix factorization based or not, default is FALSE
#' @keywords normalize
#'
#' @examples
#' dataset = list(
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50)
#' )
#' normalizeData(dataset)
#'
#' @export

## NEED TO MODIFY FOR MATRIX TYPES
normalizeData <- function(dataset, enable_normalization = TRUE, column_sum_normalization = TRUE, nonnegative_normalization = FALSE){
    if(enable_normalization){
        if(nonnegative_normalization){
          if(is.list(dataset)){
            dataset = lapply(dataset, FUN = function(x){
              cs=colSums(x);cs=cs/min(cs);sweep(x,2,cs,"/")
            })
          }else{
            cs=colSums(dataset);cs=cs/min(cs);sweep(dataset,2,cs,"/")
          }
          return(dataset)
        }
        if(column_sum_normalization){
          if(is.list(dataset)){
            dataset = lapply(dataset, FUN = function(x){
                cs=colSums(x);cs=cs/min(cs)
                t(scale(t(sweep(x,2,cs,"/")),scale = FALSE))
            })
          }else{
            cs=colSums(dataset);cs=cs/min(cs)
            t(scale(t(sweep(dataset,2,cs,"/")),scale = FALSE))
          }
        }else{
          if(is.list(dataset)){
            dataset = lapply(dataset, FUN = function(x){
              scale(t(scale(t(x), scale = FALSE)))
            })
          } else{
            dataset = scale(t(scale(t(dataset), scale = FALSE)))
          }
        }
    }else{
        dataset = lapply(dataset, FUN = function(x){
            t(scale(t(x), scale = FALSE))
        })
    }
    return(dataset)
}

#' Min-Max Normalization Function
#'
#' This function performs min-max normalization on a vector.
#' 
#' @param x A numeric vector to be normalized.
#' @return A normalized version of the input vector with values scaled between 0 and 1.
#' @examples
#' 
#' list_score = matrix(runif(5000, 1, 2), nrow = 100, ncol = 50)
#  list_score <- t(apply(t(list_score), 2, min_max_normalization))
#'
#' @export

min_max_normalization <- function(x) {
  if (max(x) - min(x) == 0) {
    return(x)  # Or return(rep(0, length(x))) if normalization to 0 is desired
  }
  (x - min(x)) / (max(x) - min(x))
}
