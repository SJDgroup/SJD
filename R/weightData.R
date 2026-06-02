#' Weighting Data Set
#'
#' To weight each data set based on input weighting vector
#'
#' @param dataset A list of data sets
#' @param weighting A vector of weighting constant for each data set
#'
#' @return A list of weighted data sets
#'
#' @keywords weighting
#'
#' @examples
#' dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' weighting = c(1, 2)
#' weighted_dataset = weightData(dataset, weighting)
#'
#' @export

weightData <- function(dataset, weighting){
    if(is.null(weighting)){
        return(dataset)
    }

    for(i in 1 : length(dataset)){
        dataset[[i]] = dataset[[i]] * weighting[i]
    }

    return(dataset)
}
