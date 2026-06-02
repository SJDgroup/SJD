#' Two-staged Independent LCA and automatic rank selection
#'
#' Two-staged decomposition of several matrices with Independent LCA,
#' twoStageLCA is first performed on the data,
#' the rank selection procedure is automatic based on BEMA.
#' Then, fastICA is implemented on the score to extract the independent
#' components.
#'
#' @param dataset A list of dataset to be analyzed
#' @param group A list of grouping of the datasets, indicating the relationship between datasets
#' @param weighting Weighting of each dataset, initialized to be NULL
#' @param threshold The threshold used to cutoff the eigenvalues
#' @param backup A backup variable, which permits the overselection of the components by BEMA
#' @param total_number Total number of components will be extracted, if default value is set to NA, then BEMA will be used.
#' @param plotting A boolean value to determine whether to plot the scree plot or not, default to be False
#' @param proj_dataset The dataset(s) to be projected on. 
#' @param proj_group A listed of boolean combinations indicating which groupings should be used for each projected dataset.The length of proj_group should match the length of proj_dataset, and the length of each concatenated boolean combination should match the length of the parameter group.
#' @param enable_normalization An argument to decide whether to use normalizaiton or not,  default is TRUE
#' @param column_sum_normalization An argument to decide whether to use column sum normalization or not, default it FALSE
#' @param screen_prob A vector of probabilies for genes to be chosen
#'
#' @importFrom RSpectra svds
#' @importFrom fastICA fastICA
#'
#' @return A list contains the component and the score of each dataset on every component after twoStageiLCA.rank algorithm
#'
#' @keywords two-staged, rank, indepentdent, LCA
#'
#' @examples
#' dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' group = list(c(1, 2, 3, 4), c(1, 2), c(3, 4), c(1, 3), c(2, 4), c(1), c(2), c(3), c(4))
#' threshold = c(3, 1.5, 1.5, 1.5, 1.5, 0.5, 0.5, 0.5, 0.5)
#' res_twoStageiLCA.rank = twoStageiLCA.rank(
#' dataset,
#' group,
#' threshold = threshold)
#'
#' @export

twoStageiLCA.rank <- function(dataset, group, weighting = NULL, total_number = NULL, threshold, backup = 0, plotting = FALSE, proj_dataset = NULL, proj_group = NULL, enable_normalization = TRUE, column_sum_normalization = FALSE, screen_prob = NULL){

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

    ## Parameters to be initialized
    N = length(dataset)
    K = length(group)
    p = nrow(dataset[[1]])
    N_dataset = unlist(lapply(dataset, ncol))

    ## Determine the number of components

    ## compute the component for each dataset
    data_comp_total = list()

    ## compute the total number of component each dataset has
    data_comp_num = rep(0, N)

    ## Using BEMA to extract the number of components
    for(i in 1 : N){
        svd_temp = svd(dataset[[i]])
        if(is.null(total_number)){
            data_comp_num[i] = BEMA(svd_temp$d^2 / ncol(dataset[[i]]), p = nrow(dataset[[i]]), n = ncol(dataset[[i]])) + backup
        }else{
            data_comp_num[i] = total_number[i] + backup
        }
        data_comp_total[[i]] = svd_temp$u[, 1 : data_comp_num[i]]
    }

    ## Output the component and scores
    list_component = list()
    list_score = list()
    for(j in 1 : N){
        list_score[[j]] = list()
    }

    for(i in 1 : K){
        list_component[[i]] = matrix(0, nrow = p, ncol = 1)
        for(j in 1 : N){
            list_score[[j]][[i]] = matrix(0, nrow = 1, ncol = N_dataset[j])
        }
    }

    ## compute the components sequentially
    for(i in 1 : K){
        temp_comp = c()
        for(j in group[[i]]){
            temp_comp = cbind(temp_comp, data_comp_total[[j]])
        }

        ## Orthogonalize the component
        if(i >= 2){
            for(j in 1 : (i - 1)){
                temp_comp = temp_comp - list_component[[j]] %*% (t(list_component[[j]]) %*% temp_comp)
            }
        }

        temp_comp_svd = svd(temp_comp)
        if(plotting){
            plot(temp_comp_svd$d^2, xlab = "index of eigenvalue", ylab = "eigenvalue", main = paste0("Group of Dataset: ", toString(group[[i]])), type = "o")
        }

        index = which(temp_comp_svd$d^2 > threshold[i])
        if(length(index) > 0){
            list_component[[i]] = matrix(temp_comp_svd$u[, 1 : length(index)], nrow = p)
            for(j in 1 : N){
                data_comp_total[[j]] = data_comp_total[[j]] - list_component[[i]] %*% (t(list_component[[i]]) %*% data_comp_total[[j]])
                # data_comp_total[[j]] = svd(data_comp_total[[j]])$u[, 1 : (ncol(data_comp_total[[j]]) - comp_num[i])]
            }
        }
    }

    ## compute the score for each dataset
    for(i in 1 : K){
        for(j in 1 : N){
            if (j %in% group[[i]]){
                list_score[[j]][[i]] = t(list_component[[i]]) %*% dataset[[j]]
            }
        }
    }

    ## Conduct ICA on the extracted Scores
    for(i in 1 : K){
        for(j in 1 : N){
            if(j %in% group[[i]] & nrow(list_score[[j]][[i]]) >= 2){
                temp_ica_res = fastICA(t(list_score[[j]][[i]]), n.comp = nrow(list_score[[j]][[i]]))
                list_score[[j]][[i]] = t(temp_ica_res$S)
            }
        }
    }

    ## Assign name for components
    list_component = compNameAssign(list_component, group_name)
    list_component = geneNameAssign(list_component, gene_name)
    list_score = scoreNameAssign(list_score, dataset_name, group_name)
    list_score = sampleNameAssign(list_score, sample_name)
    list_score = filterNAValue(list_score, dataset, group)

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
