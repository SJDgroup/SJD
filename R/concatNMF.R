#' Concatenated decomposition with Nonnegative Matrix Factorization
#'
#' Concatenated decomposition of several matrices with Nonnegative Matrix Factorization (NMF)
#'
#' @param dataset A list of dataset to be analyzed
#' @param group A list of grouping of the datasets, indicating the relationship between datasets
#' @param comp_num A vector indicates the dimension of each compoent
#' @param weighting Weighting of each dataset, initialized to be NULL
#' @param perturbation the perturbation of the 0 element in the analysis
#' @param proj_dataset The dataset(s) to be projected on. 
#' @param proj_group A listed of boolean combinations indicating which groupings should be used for each projected dataset.The length of proj_group should match the length of proj_dataset, and the length of each concatenated boolean combination should match the length of the parameter group.
#' @param enable_normalization An argument to decide whether to use normalizaiton or not,  default is TRUE
#' @param column_sum_normalization An argument to decide whether to use column sum normalization or not, default it FALSE
#' @param screen_prob A vector of probabilies for genes to be chosen
#'
#' @importFrom NMF nmf
#'
#' @return A list contains the component and the score of each dataset on every component after concatNMF algorithm
#'
#' @keywords pairwise, NMF
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
#' res_concatNMF = concatNMF(
#' dataset,
#' group,
#' comp_num,
#' proj_dataset = proj_dataset,
#' proj_group = proj_group)
#'
#'
#' @export

concatNMF = function(dataset, group, comp_num, weighting = NULL, perturbation = 0.0001, proj_dataset = NULL, proj_group = NULL, enable_normalization = TRUE, column_sum_normalization = FALSE, screen_prob = NULL,    max_ite = 1000,
                      max_err = 0.0001){

    ## Obtain names for dataset, gene and samples
    dataset_name = datasetNameExtractor(dataset)
    gene_name = geneNameExtractor(dataset)
    sample_name = sampleNameExtractor(dataset)
    group_name = groupNameExtractor(group)

    dataset = frameToMatrix(dataset)
    if(!is.null(screen_prob)){
        dataset = geneScreen(dataset, screen_prob)
    }

    dataset = normalizeData(dataset, enable_normalization, column_sum_normalization, nonnegative_normalization = TRUE)
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

    ## Extract pairwise NMF from the datasets
    for(i in 1 : K){
        temp_dat = c()
        temp_sample_n = c()
        for(j in group[[i]]){
            temp_dat = cbind(temp_dat, dataset[[j]])
            temp_sample_n = c(temp_sample_n, ncol(dataset[[j]]))
        }
        nmf_temp = nmf(temp_dat + perturbation, comp_num[i])
        list_component[[i]] = nmf_temp@fit@W
        for(j in 1 : length(group[[i]])){
            list_score[[group[[i]][j]]][[i]] = nmf_temp@fit@H[, ifelse(j == 1, 1, sum(temp_sample_n[1 : (j - 1)]) + 1) : sum(temp_sample_n[1 : j])]
        }
    }

    ## Assign name for components
    list_component = compNameAssign(list_component, group_name)
    list_component = geneNameAssign(list_component, gene_name)
    list_score = scoreNameAssign(list_score, dataset_name, group_name)
    list_score = sampleNameAssign(list_score, sample_name)
    list_score = filterNAValue(list_score, dataset, group)
    list_score = rebalanceData(list_score, group, dataset)
    
    
    ## Max-normalize sample scores while preserving reconstruction
    eps = 1e-12
    
    for (i in seq_len(K)) {
      for (k in seq_len(comp_num[i])) {
        
        # Find the maximum score across all datasets
        # participating in this component group
        component_max = 0
        
        for (j in group[[i]]) {
          component_max = max(
            component_max,
            max(list_score[[j]][[i]][k, ], na.rm = TRUE)
          )
        }
        
        if (component_max > eps) {
          
          # Normalize sample scores
          for (j in group[[i]]) {
            list_score[[j]][[i]][k, ] =
              list_score[[j]][[i]][k, ] / component_max
          }
          
          # Compensate by scaling gene embeddings
          list_component[[i]][, k] =
            list_component[[i]][, k] * component_max
        }
      }
    }


    ## Project score
    if(!is.null(proj_dataset)){
      if (length(proj_group) != length(list_component)){
        stop("Error:length of proj_group should equal to length of list_component.")
      }
      
      proj_dataset = as.matrix(proj_dataset)
      
      proj_component = list_component
      
      if(!is.null(rownames(proj_dataset))){
        CMNgenes=rownames(proj_dataset)[rownames(proj_dataset)%in%rownames(proj_component[[1]])]
        common_genes.DATA = match(CMNgenes,rownames(proj_dataset))
        common_genes.COMP = match(CMNgenes,rownames(proj_component[[1]]))
        orig_gene_length = nrow(proj_dataset)
        proj_dataset = proj_dataset[common_genes.DATA,]
        for (i in 1:length(proj_component)) {
          proj_component[[i]] = proj_component[[i]][common_genes.COMP,]
        }
        print(paste("Input", orig_gene_length, "genes in proj_dataset, found", nrow(proj_dataset), "genes in common."))
      }
      
      
      p = nrow(proj_dataset)
      comp_num = unlist(lapply(proj_component, function(x) if (length(x) > 0) ncol(x)))
      
      # Initialize the matrix W with zeros
      W = matrix(0, nrow = p, ncol = sum(comp_num))
      K = length(proj_group)
      # Reconstruct matrix W using proj_component and comp_num
      start_col = 1
      for (i in 1:K) {
        end_col = start_col + comp_num[i] - 1
        W[, start_col:end_col] = proj_component[[i]]
        start_col = end_col + 1
      }
      
      
      proj_sample_name = sampleNameExtractor(proj_dataset)
      group_name = names(proj_component)
      proj_dataset = normalizeData(proj_dataset, enable_normalization, column_sum_normalization, nonnegative_normalization = TRUE)
      
      M = sum(comp_num)
      col = ncol(proj_dataset)
      
      
      ## Initialize the W and H for Nonnegative Matrix Factorization
      max_element = -Inf
      min_element = Inf
      max_element = max(max_element, max(proj_dataset))
      min_element = min(min_element, min(proj_dataset))
      
      
      ## Initialize the values of W and H
      X = proj_dataset
      
      
      H = NULL
      
      for (i in seq_len(K)) {
        
        if (proj_group[i]) {
          H_temp = matrix(
            runif(col * comp_num[i], min_element, max_element),
            nrow = comp_num[i],
            ncol = col
          )
        } else {
          H_temp = matrix(
            0,
            nrow = comp_num[i],
            ncol = col
          )
        }
        
        H = rbind(H, H_temp)
      }
      
      
      ## Iteratively estimate the NMF with Euclidean distance
      error_out = c()
      
      for(ite in 1 : max_ite){
        ## H is score (sample matrix S)
        H = H * (t(W) %*% X) /
          pmax(t(W) %*% W %*% H, 1e-12) # W is proj_component (common gene matrix G)
        
        H[which(is.na(H))] = 0
        
        error_out = c(error_out, sum((X - W %*% H)^2))
        ## Break when the error difference is small
        if(length(error_out) >= 2 && abs(error_out[length(error_out)] - error_out[length(error_out) - 1]) / abs(error_out[length(error_out) - 1]) <= max_err){
          break
        }
      }
      
      
      proj_list_score= list()
      for(j in 1 : length(proj_group)){ # j goes from 1 to 3 since we have 3 groups
        if(proj_group[j]){
          proj_list_score[[j]] = H[ifelse(j == 1, 1, cumsum(comp_num)[j - 1] + 1) : cumsum(comp_num)[j],  1:col]
        } else {
          proj_list_score[[j]] = matrix(0, nrow = comp_num[j], ncol = col)
        }
      }
      
      proj_list_score = scoreNameAssignProj(proj_list_score, group_name)
      proj_list_score = sampleNameAssignProj(proj_list_score, proj_sample_name)
      # for(i in 1 : length(proj_group)){
      #   if(proj_group[j]){
      #     proj_list_score[[i]] = proj_list_score[[i]] * sqrt(ncol(proj_dataset))
      #   }
      # }
    } else {
      proj_list_score = NULL
    }

    return(list(linked_component_list = list_component, score_list = list_score, proj_score_list = proj_list_score))
}
