#' Joint decomposition with Independent Component Analysis
#'
#' Joint decomposition of several linked matrices with Independent Component Analysis (ICA)
#'
#' @param dataset A list of dataset to be analyzed
#' @param group A list of grouping of the datasets, indicating the relationship between datasets
#' @param comp_num A vector indicates the dimension of each compoent
#' @param weighting Weighting of each dataset, initialized to be NULL
#' @param max_ite The maximum number of iterations for the jointPCA algorithms to run, default value is set to 100
#' @param max_err The maximum error of loss between two iterations, or the program will terminate and return, default value is set to be 0.0001
#' @param proj_dataset The dataset(s) to be projected on. 
#' @param proj_group A listed of boolean combinations indicating which groupings should be used for each projected dataset.The length of proj_group should match the length of proj_dataset, and the length of each concatenated boolean combination should match the length of the parameter group.
#' @param enable_normalization An argument to decide whether to use normalizaiton or not,  default is TRUE
#' @param column_sum_normalization An argument to decide whether to use column sum normalization or not, default it FALSE
#' @param screen_prob A vector of probabilies for genes to be chosen
#'
#' @importFrom fastICA fastICA
#'
#' @return A list contains the component and the score of each dataset on every component after jointPCA algorithm
#'
#' @keywords Joint, ICA
#'
#' @examples
#' dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' group = list(c(1,2,3,4), c(1,2), c(3,4), c(1,3), c(2,4), c(1), c(2), c(3), c(4))
#' comp_num = c(2,2,2,2,2,2,2,2,2)
#' proj_dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' proj_group = list(c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))
#' res_jointICA = jointICA(
#' dataset,
#' group,
#' comp_num,
#' proj_dataset = proj_dataset,
#' proj_group = proj_group)
#'
#'
#' @export

jointICA <- function(dataset, group, comp_num, weighting = NULL, max_ite = 100, max_err = 0.0001, proj_dataset = NULL, proj_group = NULL, enable_normalization = TRUE, column_sum_normalization = FALSE, screen_prob = NULL){
    jointPCA_res = jointPCA(dataset, group, comp_num, weighting, max_ite, max_err)

    ## Obtain names for dataset, gene and samples
    dataset_name = datasetNameExtractor(dataset)
    gene_name = geneNameExtractor(dataset)
    sample_name = sampleNameExtractor(dataset)
    group_name = groupNameExtractor(group)

    dataset = frameToMatrix(dataset)
    if(!is.null(screen_prob)){
        dataset = geneScreen(dataset, screen_prob)
    }

    dataset = normalizeData(dataset, enable_normalization, column_sum_normalization)
    dataset = balanceData(dataset)
    dataset = weightData(dataset, weighting)

    ## Parameters to be initialized
    N = length(dataset)
    K = length(group)
    M = sum(comp_num)
    p = nrow(dataset[[1]])
    N_dataset = unlist(lapply(dataset, ncol))

    ## Output the component and scores
    list_component = list()
    list_score = list()
    for(j in 1 : N){
        list_score[[j]] = list()
    }

    for(i in 1 : K){
        list_component[[i]] = matrix(0, nrow = p, ncol = comp_num[i])
        for(j in 1 : N){
            list_score[[j]][[i]] = matrix(0, nrow = comp_num[i], ncol = N_dataset[j])
        }
    }

    ## Applying joint ICA
    for(i in 1 : K){
        dat_temp = c()
        n_sample_temp = c()
        for(j in group[[i]]){
            dat_temp = cbind(dat_temp, jointPCA_res$score_list[[j]][[i]] / sqrt(N_dataset[j]))
            n_sample_temp = c(n_sample_temp, N_dataset[j])
        }
        ica_temp = fastICA(t(dat_temp), comp_num[i])
        list_component[[i]] = jointPCA_res$linked_component_list[[i]] %*% t(ica_temp$A)
        for(j in 1 : length(group[[i]])){
            list_score[[group[[i]][j]]][[i]] = t(ica_temp$S)[, ifelse(j == 1, 1, sum(n_sample_temp[1 : (j - 1)]) + 1) : sum(n_sample_temp[1 : j])]
        }
    }

    ## Assign name for components
    list_component = compNameAssign(list_component, group_name)
    list_component = geneNameAssign(list_component, gene_name)
    list_score = scoreNameAssign(list_score, dataset_name, group_name)
    list_score = sampleNameAssign(list_score, sample_name)
    list_score = filterNAValue(list_score, dataset, group)
    list_score = rebalanceData(list_score, group, dataset)
    list_score = pveMultiple(dataset, group, comp_num, list_score, list_component)

    ## Project score
    proj_list_score = list()
    if(!is.null(proj_dataset)){
        proj_sample_name = sampleNameExtractor(proj_dataset)
        proj_dataset_name = datasetNameExtractor(proj_dataset)

        proj_dataset = frameToMatrix(proj_dataset)
        proj_dataset = normalizeData(proj_dataset, enable_normalization, column_sum_normalization)
        proj_dataset = balanceData(proj_dataset)

        for(i in 1 : length(proj_dataset)){
            proj_list_score[[i]] = list()
            for(j in 1 : length(proj_group[[i]])){
                if(proj_group[[i]][j]){
                    proj_list_score[[i]][[j]] =  t(list_component[[j]]) %*% proj_dataset[[i]]
                }else{
                    proj_list_score[[i]][[j]] = NA
                }
            }
        }


        proj_list_score = scoreNameAssign(proj_list_score, proj_dataset_name, group_name)
        proj_list_score = sampleNameAssign(proj_list_score, proj_sample_name)
    }


    return(list(linked_component_list = list_component, score_list = list_score, proj_score_list = proj_list_score))
}
