#' Single Data Set Decomposition with Independent Component Analysis
#'
#' Apply ICA (Independent Component Analysis) to a single data set
#'
#' @param dataset A dataframe/matrix to be decomposed
#' @param comp_num Number of ICs to be extracted
#' @param weighting Weighting of each dataset, initialized to be NULL
#' @param enable_normalization An argument to decide whether to use normalizaiton or not,  default is TRUE
#' @param column_sum_normalization An argument to decide whether to use column sum normalization or not, default it FALSE
#' @param screen_prob A vector of probabilies for genes to be chosen
#'
#' @importFrom fastICA fastICA
#'
#' @return A list of scores and component
#'
#' @keywords separate analysis, ICA
#'
#' @examples
#' dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' comp_num = 2
#' res_sepICA = sepICA(dataset, comp_num)
#'
#' @export

sepICA <- function(dataset, comp_num, weighting = NULL, enable_normalization = TRUE, column_sum_normalization = FALSE, screen_prob = NULL){
    sepPCA_res = sepPCA(dataset, comp_num)

    ## Obtain names for dataset, gene and samples
    dataset_name = datasetNameExtractor(dataset)
    gene_name = geneNameExtractor(dataset)
    sample_name = sampleNameExtractor(dataset)

    dataset = frameToMatrix(dataset)

    if(!is.null(screen_prob)){
        dataset = geneScreen(dataset, screen_prob)
    }

    dataset = normalizeData(dataset, enable_normalization, column_sum_normalization)
    dataset = weightData(dataset, weighting)

    N = length(dataset)

    list_component = list()
    list_score = list()

    for(i in 1 : N){
        ica_temp = fastICA(t(sepPCA_res$score_list[[i]]), comp_num[i])
        component = sepPCA_res$linked_component_list[[i]] %*% t(ica_temp$A)
        score = t(ica_temp$S)
        list_component[[i]] = component
        list_score[[i]] = score
    }

    ## Assign name for components
    list_component = compNameAssignSep(list_component, dataset_name)
    list_component = geneNameAssign(list_component, gene_name)
    list_score = scoreNameAssignSep(list_score, dataset_name)
    list_score = sampleNameAssignSep(list_score, sample_name)
    list_score = pveSep(dataset, list_score, list_component)

    return(list(linked_component_list = list_component, score_list = list_score))
}
