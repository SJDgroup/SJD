#' Transform data.frame to matrix
#'
#' Function to transform a data frame into a matrix
#'
#' @param dataset A list of data sets
#'
#' @return A list of matrix transformed from the list of data frames
#'
#' @keywords dataframe, matrix
#'
#' @examples
#' dataset = list(
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50)
#' )
#' frameToMatrix(dataset)
#'
#' @export


frameToMatrix <- function(dataset){
    dataset = lapply(dataset, FUN = function(x){as.matrix(x)})
    return(dataset)
}
