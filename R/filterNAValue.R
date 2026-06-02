#' Filter NA Value
#'
#' Assign NA values to scores not in the group for each data set
#'
#' @param list_score A list of scores extracted in the analysis
#' @param dataset A list of datasets used in the analysis
#' @param group A list of group assignments
#'
#' @return A list of scores by assigning NA to scores not in groups
#'
#' @keywords NA, filter
#'
#' @examples
#'
#' x = list(matrix(c(1 : 4), nrow = 2), matrix(c(1 : 4), nrow = 2))
#' y = list(matrix(c(1 : 4), nrow = 2), matrix(c(1 : 4), nrow = 2))
#' list_score = list(x, y)
#' dataset = c(1, 2)
#' group = list(c(1), c(2))
#' filterNAValue(list_score, dataset, group)
#'
#' @export

filterNAValue <- function(list_score, dataset, group){
    for(i in 1 : length(dataset)){
        for(j in 1 : length(group)){
            if(!(i %in% group[[j]])){
                list_score[[i]][[j]] = NA
            }
        }
    }
    return(list_score)
}
