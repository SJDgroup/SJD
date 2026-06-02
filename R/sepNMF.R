#' Single Data Set Decomposition with Nonnegative Matrix Factorization
#'
#' Apply NMF (Nonnegative Matrix Factorization) to a single data set
#'
#' @param dataset A dataframe/matrix to be decomposed
#' @param comp_num Number of NMFs to be extracted
#' @param weighting Weighting of each dataset, initialized to be NULL
#' @param perturbation A small perturbation to ensure nmf works well
#' @param enable_normalization An argument to decide whether to use normalizaiton or not,  default is TRUE
#' @param column_sum_normalization An argument to decide whether to use column sum normalization or not, default it FALSE
#' @param screen_prob A vector of probabilies for genes to be chosen
#'
#' @importFrom NMF nmf
#'
#' @return A list of scores and component
#'
#' @keywords separate analysis, NMF
#'
#' @examples
#' dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' comp_num = 2
#' res_sepNMF = sepNMF(dataset, comp_num)
#'
#' @export

sepNMF <- function(dataset, comp_num, weighting = NULL, perturbation = 0.0001, enable_normalization = TRUE, column_sum_normalization = FALSE, screen_prob = NULL){

    ## Obtain names for dataset, gene and samples
    dataset_name = datasetNameExtractor(dataset)
    gene_name = geneNameExtractor(dataset)
    sample_name = sampleNameExtractor(dataset)

    dataset = frameToMatrix(dataset)

    if(!is.null(screen_prob)){
        dataset = geneScreen(dataset, screen_prob)
    }

    dataset = normalizeData(dataset, enable_normalization, column_sum_normalization, nonnegative_normalization = TRUE)
    dataset = weightData(dataset, weighting)

    N = length(dataset)

    list_component = list()
    list_score = list()

    for(i in 1 : N){
        nmf_temp = nmf(dataset[[i]] + perturbation, comp_num[i])
        component = nmf_temp@fit@W
        score = nmf_temp@fit@H
        list_component[[i]] = component
        list_score[[i]] = score
    }

    ## Assign name for components
    list_component = compNameAssignSep(list_component, dataset_name)
    list_component = geneNameAssign(list_component, gene_name)
    list_score = scoreNameAssignSep(list_score, dataset_name)
    list_score = sampleNameAssignSep(list_score, sample_name)

    return(list(linked_component_list = list_component, score_list = list_score))
}
