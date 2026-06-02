#' Gene Screen
#'
#' Screen genes based on their standard deviation
#'
#' @param dataset A list of dataset to be analyzed
#' @param screen_prob A vector of probabilies for genes to be chosen
#'
#' @return A list contains the component and the score of each dataset on every component after jointPCA algorithm
#'
#' @keywords screen, gene
#'
#' @examples
#' dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#'
#' screen_prob = c(0.2, 0.2, 0.2, 0.2)
#'
#' screened_dataset = geneScreen(dataset, screen_prob)
#'
#' @export

geneScreen <- function(dataset, screen_prob){
    p = nrow(dataset[[1]])
    N = length(dataset)
    gene_index = list()

    for(i in 1 : N){
        geneSd = apply(dataset[[i]], 1, sd)
        cutoff = quantile(geneSd, screen_prob[i])
        gene_index[[i]] = which(geneSd > cutoff)
    }

    geneSelected = c(1 : p)
    for(i in 1 : N){
        geneSelected = intersect(geneSelected, gene_index[[i]])
    }

    output_dataset = list()

    for(i in 1 : N){
        output_dataset[[i]] = dataset[[i]][geneSelected, ]
    }

    return(output_dataset)
}
