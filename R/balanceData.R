#' Balance Data
#'
#' Balance Data for Concatenate and Joint Analysis
#'
#' @param dataset A list of data sets to be used
#'
#' @return A list of rebalanced data sets
#'
#' @keywords balance, dataset
#'
#' @examples
#' dataset = list(matrix(c(1 : 8), nrow = 2), matrix(1 : 6, nrow = 2))
#' balanceData(dataset)
#'
#' @export

balanceData <- function(dataset){
    for(i in 1 : length(dataset)){
        dataset[[i]] = dataset[[i]] / sqrt(ncol(dataset[[i]]))
    }
    return(dataset)
}



#' Rebalance Data
#'
#' Rebalance scores based on the balanced data set
#'
#' @param list_score A list of scores extracted in the analysis
#' @param group A list of group assignment
#' @param dataset A list of dataset in analysis
#'
#' @return A list of rebalanced scores
#'
#' @examples
#' x = list(matrix(c(1 : 4), nrow = 2), matrix(c(1 : 4), nrow = 2))
#' y = list(matrix(c(1 : 4), nrow = 2), matrix(c(1 : 4), nrow = 2))
#' list_score = list(x, y)
#' group = list(c(1), c(2))
#' dataset = list(matrix(c(1 : 8), nrow = 2), matrix(1 : 6, nrow = 2))
#' rebalanceData(list_score, group, dataset)
#'
#' @export

rebalanceData <- function(list_score, group, dataset){
    for(i in 1 : length(group)){
        for(j in 1 : length(dataset)){
            if(j %in% group[[i]]){
                list_score[[j]][[i]] = list_score[[j]][[i]] *  sqrt(ncol(dataset[[j]]))
            }
        }
    }
    return(list_score)
}
