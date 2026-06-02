#' Two-staged Independent Linked Component Analysis
#'
#' Two-staged Independent Linked Component Analysis, a generalization based on the Two-staged Independent Linked Component Analysis
#'
#' @param dataset A list of dataset to be analyzed
#' @param group A list of grouping of the datasets, indicating the relationship between datasets
#' @param comp_num A vector indicates the dimension of each compoent
#' @param weighting Weighting of each dataset, initialized to be NULL
#' @param backup A positive scalar to determine how many ICs to over select
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
#' @return A list contains the component and the score of each dataset on every component after 2siLCA algorithm
#'
#' @keywords two-staged, independent LCA
#'
#' @examples
#' dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' group = list(c(1,2,3,4), c(1,2), c(3,4), c(1,3), c(2,4), c(1), c(2), c(3), c(4))
#' comp_num = c(2,2,2,2,2,2,2,2,2)
#' proj_dataset = matrix(runif(5000, 1, 2), nrow = 100, ncol = 50)
#' proj_group = c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)
#' res_twoStageiLCA = twoStageiLCA(
#' dataset,
#' group,
#' comp_num,
#' proj_dataset = proj_dataset,
#' proj_group = proj_group)
#'
#'
#' @export

twoStageiLCA <- function(dataset, group, comp_num, weighting = NULL, backup = 0, plotting = FALSE, proj_dataset = NULL, proj_group = NULL, enable_normalization = TRUE, column_sum_normalization = FALSE, screen_prob = NULL){
    twoStageLCA_out = twoStageLCA(dataset, group, comp_num, weighting, backup)

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


    ## Conduct ICA on the extracted Scores
    ica_score = list()
    for(i in 1 : K){
        list_component[[i]] = twoStageLCA_out$linked_component_list[[i]]
        score_concat  = c()
        for(j in 1 : N){
          if(j %in% group[[i]] & nrow(list_score[[j]][[i]]) >= 2){
            score_concat = cbind(score_concat, twoStageLCA_out$score_list[[j]][[i]])
            # print(dim(twoStageLCA_out$score_list[[j]][[i]]))
            # ica_temp = fastICA(t(twoStageLCA_out$score_list[[j]][[i]]), n.comp = nrow(twoStageLCA_out$score_list[[j]][[i]]))
            # list_score[[j]][[i]] = t(ica_temp$S)
          }
        }
        # SCORES ENTER ICA?? WHY DO WE WANT TO DECOMPOSE THE SCORES from LCA?
        # The fundamental goal of ICA is to uncover hidden source signals that are mixed in the observed data.
        ica_temp = fastICA(t(score_concat), n.comp = nrow(score_concat))
        ica_score[[i]] = ica_temp
        start_index = 0
        
        for(j in 1 : N){
            if(j %in% group[[i]] & nrow(list_score[[j]][[i]]) >= 2){
              
              # NEW SCORES (SAMPLE SCORES) ARE ESTIMATED SOURCE MATRICES X = AS
              list_score[[j]][[i]] = t(ica_temp$S[(start_index + 1) : (start_index + N_dataset[j]), ]) #slice ica result
              start_index = start_index + N_dataset[j]
            }
        }
    }
        

    ## Assign name for components
    list_component = compNameAssign(list_component, group_name)
    list_component = geneNameAssign(list_component, gene_name)
    list_score = scoreNameAssign(list_score, dataset_name, group_name)
    list_score = sampleNameAssign(list_score, sample_name)
    list_score = filterNAValue(list_score, dataset, group)
    list_score = pveMultiple(dataset, group, comp_num, list_score, list_component)

    ## Project score
    proj_list_score = list()
    if(!is.null(proj_dataset)){
        proj_sample_name = sampleNameExtractor(proj_dataset)
        group_name = names(list_component) 
        proj_dataset = normalizeData(proj_dataset, enable_normalization, column_sum_normalization)

        for(j in 1 : length(proj_group)){
            if(proj_group[[j]]){
                proj_list_score[[j]] = t(t((t(list_component[[j]]) %*% proj_dataset)) %*% ica_score[[j]]$K %*% ica_score[[j]]$W)
            }else{
                proj_list_score[[j]] = NA
            }
        }

        proj_list_score = scoreNameAssignProj(proj_list_score, group_name)
        proj_list_score = sampleNameAssignProj(proj_list_score, proj_sample_name)
    }


    return(list(linked_component_list = list_component, score_list = list_score, proj_score_list = proj_list_score, ica_score = ica_score))
}
