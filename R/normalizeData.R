#' Data Normalization
#'
#' Normalize data to have mean zero and std 1
#'
#' @param dataset The input list of data sets matrix
#' @param enable_normalization An argument to decide whether to use normalizaiton or not,  default is TRUE
#' @param column_sum_normalization An argument to decide whether to use column sum normalization or not, default it TRUE
#' @param nonnegative_normalization An argument to decide wehther it is nonnegative matrix factorization based or not, default is FALSE. If TRUE, each column is divided by its sum relative to the smallest column sum (data stay nonnegative); if TRUE and enable_normalization is FALSE, the data are returned unchanged.
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
    if(nonnegative_normalization){
        ## NMF-based methods: rescale columns to a common sum, which keeps the
        ## data nonnegative. With enable_normalization = FALSE the data are used
        ## as provided (e.g. already normalized by the caller).
        if(!enable_normalization){
            return(dataset)
        }
        if(is.list(dataset)){
            for(i in seq_along(dataset)){
                dataset[[i]] = columnSumScale(dataset[[i]])
            }
        }else{
            dataset = columnSumScale(dataset)
        }
        return(dataset)
    }
    if(enable_normalization){
        if(column_sum_normalization){
          if(is.list(dataset)){
            dataset = lapply(dataset, FUN = function(x){
                cs=colSums(x);cs=cs/min(cs)
                t(scale(t(sweep(x,2,cs,"/")),scale = FALSE))
            })
          }else{
              cs=colSums(dataset);cs=cs/min(cs)
              dataset = t(scale(t(sweep(dataset,2,cs,"/")),scale = FALSE))
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
        if(is.list(dataset)){
            dataset = lapply(dataset, FUN = function(x){
                t(scale(t(x), scale = FALSE))
            })
        }else{
            dataset = t(scale(t(dataset), scale = FALSE))
        }
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

#' Column-sum scaling without large temporaries
#'
#' Divides every column of \code{x} by its sum relative to the smallest column
#' sum, i.e. the same result as
#' \code{cs <- colSums(x); sweep(x, 2, cs / min(cs), "/")}, but without the
#' full-size array that \code{sweep()} builds. Columns are processed in blocks
#' of about \code{block_bytes} so only one copy of \code{x} is made.
#'
#' @param x A numeric matrix.
#' @param block_bytes Approximate size of each processed block, in bytes.
#'
#' @return The scaled matrix (same dimensions and dimnames as \code{x}).
#'
#' @keywords internal
#' @noRd

columnSumScale <- function(x, block_bytes = 64 * 1024^2){
    if(is.data.frame(x)){
        x <- as.matrix(x)
    }
    p <- nrow(x)
    n <- ncol(x)
    if(n == 0L){
        return(x)
    }
    cs <- colSums(x)
    cs <- cs / min(cs)
    block <- max(1L, as.integer(block_bytes / (8 * max(p, 1L))))
    for(start in seq(1L, n, by = block)){
        idx <- start:min(start + block - 1L, n)
        x[, idx] <- x[, idx, drop = FALSE] / rep(cs[idx], each = p)
        ## Free the block temporaries now (cheap minor collection); otherwise
        ## they pile up as garbage because R schedules collections relative to
        ## the (large) heap size.
        gc(verbose = FALSE, full = FALSE)
    }
    x
}
