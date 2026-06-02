#' Joint Decomposition with Principal Component Analysis
#'
#' Joint decomposition of several linked matrices with Principal Component Analysis (PCA)
#'
#' @param dataset A list of dataset to be analyzed
#' @param group A list of grouping of the datasets, indicating the relationship between datasets
#' @param comp_num A vector indicates the dimension of each compoent
#' @param weighting Weighting of each dataset, initialized to be NULL
#' @param max_ite The maximum number of iterations for the jointPCA algorithms to run, default value is set to 100
#' @param max_err The maximum error of loss between two iterations, or the program will terminate and return, default value is set to be 0.001
#' @param proj_dataset The dataset(s) to be projected on. 
#' @param proj_group A listed of boolean combinations indicating which groupings should be used for each projected dataset.The length of proj_group should match the length of proj_dataset, and the length of each concatenated boolean combination should match the length of the parameter group.
#' @param enable_normalization An argument to decide whether to use normalizaiton or not,  default is TRUE
#' @param column_sum_normalization An argument to decide whether to use column sum normalization or not, default it FALSE
#' @param screen_prob A vector of probabilies for genes to be chosen
#'
#' @importFrom stats runif
#'
#' @return A list contains the component and the score of each dataset on every component after jointPCA algorithm
#'
#' @keywords joint, PCA
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
#' res_jointPCA = jointPCA(
#' dataset,
#' group,
#' comp_num,
#' proj_dataset = proj_dataset,
#' proj_group = proj_group)
#'
#'
#' @export

jointPCA <- function(dataset, group, comp_num, weighting = NULL, max_ite = 100, max_err = 0.0001, proj_dataset = NULL, proj_group = NULL, enable_normalization = TRUE, column_sum_normalization = FALSE, screen_prob = NULL){

    ## Obtain names for dataset, gene and samples
    dataset_name = datasetNameExtractor(dataset)
    gene_name = geneNameExtractor(dataset)
    sample_name = sampleNameExtractor(dataset)
    group_name = groupNameExtractor(group)

    ## Preprocess datasets
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

    ## Combine the dataset into a huge one
    combine_data <- c()
    for(i in 1 : N){
        combine_data = cbind(combine_data, dataset[[i]])
    }

    ## List to store the random scores and initialize the scores
    list_score = list()
    for(i in 1 : N){
        list_score[[i]] = list()
    }
    for(i in 1 : N){
        for(j in 1 : K){
            if(i %in% group[[j]]){
                list_score[[i]][[j]] = matrix(runif(comp_num[j] * N_dataset[i]), nrow = comp_num[j])
            }else{
                list_score[[i]][[j]] = matrix(0, nrow = comp_num[j], ncol = N_dataset[i])
            }
        }
    }

    ## Initialize the loss
    loss = 0

    ## Start the Alternative Projection
    for(t in 1 : max_ite){
        matrix_score = c()
        for(i in 1 : N){
            temp_score = c()
            for(j in 1 : K){
                temp_score = rbind(temp_score, list_score[[i]][[j]])
            }
            matrix_score = cbind(matrix_score, temp_score)
        }

        ## Apply Procrustes projection to obtain the linkedin component
        linked_component = Procrustes(combine_data, matrix_score)
        list_component = list()
        index = 1
        for(j in 1 : K){
            list_component[[j]] = linked_component[, index : (index + comp_num[j] - 1)]
            index = index + comp_num[j]
        }

        ## Compute the random scores for every dataset
        for(i in 1 : N){
            for(j in 1 : K){
                if(i %in% group[[j]]){
                    list_score[[i]][[j]] = t(list_component[[j]]) %*% dataset[[i]]
                }else{
                    list_score[[i]][[j]] = matrix(0, nrow = comp_num[j], ncol = N_dataset[i])
                }
            }
        }

        loss_current = sum((combine_data - linked_component %*% matrix_score)^2)
        if(abs(loss[length(loss)] - loss_current) < max_err){
            list_component = compNameAssign(list_component, group_name)
            list_component = geneNameAssign(list_component, gene_name)
            list_score = scoreNameAssign(list_score, dataset_name, group_name)
            list_score = sampleNameAssign(list_score, sample_name)

            return(list(linked_component_list = list_component, score_list = list_score))
        }
        loss = c(loss, loss_current)
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


    return(list(linked_component_list = list_component, score_list = list_score, loss = loss, proj_score_list = proj_list_score))
}

#' Procrustes Projection
#'
#' Procrustes projection function to solve the procrustes problem ||A - UB||^2 U^T U = I
#'
#' @param A The input matrix A as target
#' @param B the input matrix B as basis
#'
#' @return The procrustes matrix U
#'
#' @keywords procrustes
#'
#' @export

Procrustes <- function(A, B){
    C = A %*% t(B)
    svd.temp = svd(C)
    return(svd.temp$u %*% t(svd.temp$v))
}

