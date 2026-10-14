#' Joint Flexible Decomposition with Nonnegative Matrix Factorization
#'
#' Joint flexible decomposition of several linked matrices with Nonnegative Matrix Factorization (NMF)
#' It combines the jointNMF function and the flexNMF function, allowing flexible decomposition of several matrices with known patterns (gene embeddings and/or sample scores).
#' It is based on the MSE loss, proposed by Lee, Daniel D., and H. Sebastian Seung. "Learning the parts of objects by non-negative matrix factorization." Nature 401.6755 (1999): 788-791.
#'
#' @param dataset A list of dataset to be analyzed
#' @param group A list of grouping of the datasets, indicating the relationship between datasets
#' @param comp_num A vector indicates the dimension of each compoent
#' @param gene_score The known gene embedding columns/patterns to be fixed
#' @param sample_score The known sample score columns/patterns to be fixed
#' @param weighting Weighting of each dataset, initialized to be NULL
#' @param max_ite The maximum number of iterations for the jointNMF algorithms to run, default value is set to 100
#' @param max_err The maximum error of loss between two iterations, or the program will terminate and return, default value is set to be 0.0001
#' @param enable_normalization An argument to decide whether to use normalizaiton or not,  default is TRUE
#' @param column_sum_normalization An argument to decide whether to use column sum normalization or not, default it FALSE
#' @param screen_prob A vector of probabilies for genes to be chosen
#'
#' @return A list contains the component and the score of each dataset on every component after jointFlexNMF algorithm
#'
#' @keywords joint, flexible, NMF, decomposition
#'
#' @examples
#' dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' group = list(c(1,2,3,4), c(1,2), c(3,4), c(1,3), c(2,4), c(1), c(2), c(3), c(4))
#' comp_num = c(2,2,2,2,2,2,2,2,2)
#' gene_score = matrix(runif(200, 1, 2), nrow = 100, ncol = 2)
#' sample_score = matrix(runif(400, 1, 2), nrow = 2, ncol = 200)
#' res_jointFlexNMF = jointFlexNMF(
#' dataset,
#' group,
#' comp_num,
#' gene_score = gene_score,
#' sample_score = sample_score)
#'
#'
#' @export

jointFlexNMF <- function(dataset, group, comp_num, gene_score = NULL, sample_score = NULL, weighting = NULL, max_ite = 1000, max_err = 0.0001, enable_normalization = TRUE, column_sum_normalization = FALSE, screen_prob = NULL){

  ## Obtain names for dataset, gene and samples
  dataset_name = datasetNameExtractor(dataset)
  gene_name = geneNameExtractor(dataset)
  sample_name = sampleNameExtractor(dataset)
  group_name = groupNameExtractor(group)

  ## Preprocess Dataset
  dataset = frameToMatrix(dataset)

  if(!is.null(screen_prob)){
    dataset = geneScreen(dataset, screen_prob)
  }

  dataset = normalizeData(dataset, enable_normalization, column_sum_normalization, nonnegative_normalization = TRUE)
  dataset = balanceData(dataset)
  dataset = weightData(dataset, weighting)

  ## Initialize values for the algorithm
  N = length(dataset)
  K = length(group)
  M = sum(comp_num)
  p = nrow(dataset[[1]])
  N_dataset = unlist(lapply(dataset, ncol))
  total_N = sum(N_dataset)

  ## Initialize the W and H for Nonnegative Matrix Factorization
  max_element = -Inf
  min_element = Inf
  for(i in 1 : N){
    max_element = max(max_element, max(dataset[[i]]))
    min_element = min(min_element, min(dataset[[i]]))
  }

  ## Initialize the values of W and H
  X = c()
  for(i in 1 : N){
    X = cbind(X, dataset[[i]])
  }
  W = matrix(runif(p * M, min_element, max_element), nrow = p)
  H = c()

  for(i in 1 : K){
    H_temp = c()
    for(j in 1 : N){
      if(j %in% group[[i]]){
        H_temp = cbind(H_temp, matrix(runif(N_dataset[j] * comp_num[i], min_element, max_element), nrow = comp_num[i], ncol = N_dataset[j]))
      }else{
        H_temp = cbind(H_temp, matrix(0, nrow = comp_num[i], ncol = N_dataset[j]))
      }
    }
    H = rbind(H, H_temp)
  }
  
  ## Validate fixed sample scores against group structure
  if (!is.null(sample_score)) {
    
    if (nrow(sample_score) != M || ncol(sample_score) != total_N) {
      stop("sample_score must have M rows and total_N columns.")
    }
    
    component_offset = 0
    
    for (i in 1:K) {
      
      component_indices = (component_offset + 1):(component_offset + comp_num[i])
      sample_offset = 0
      
      for (j in 1:N) {
        
        sample_indices = (sample_offset + 1):(sample_offset + N_dataset[j])
        
        # This component group should not be active in dataset j
        if (!(j %in% group[[i]])) {
          
          if (any(sample_score[component_indices, sample_indices, drop = FALSE] != 0)) {
            stop(paste0(
              "Invalid sample_score: component group ", i,
              " is not active in dataset ", j,
              ", but sample_score contains nonzero values in this block."
            ))
          }
        }
        
        sample_offset = sample_offset + N_dataset[j]
      }
      
      component_offset = component_offset + comp_num[i]
    }
  }

  ## Handle pre-defined fixed columns in W (gene_score) and fixed rows in H (sample_score).
  ## Non-zero columns of gene_score are pinned in W; non-zero rows of sample_score are pinned in H.
  ## (samples concatenated in dataset order, matching the internal H layout).
  has_fixed_W <- FALSE
  has_fixed_H <- FALSE
  non_zero_cols_W <- integer(0)
  non_zero_rows_H <- integer(0)
  

  if(!is.null(gene_score) && !is.null(sample_score)){
    
    num_patterns_g = ncol(gene_score)
    num_patterns_s = nrow(sample_score)
    
    if(num_patterns_g != num_patterns_s){
      stop("ncol(gene_score) must match nrow(sample_socre)")
    }
      

    non_zero_cols_W <- which(colSums(gene_score) > 0)
    if (length(non_zero_cols_W) > 0) {
      W[, non_zero_cols_W] <- gene_score[, non_zero_cols_W, drop = FALSE]
      has_fixed_W <- TRUE
    }
  

    non_zero_rows_H <- which(rowSums(sample_score) > 0)
    if (length(non_zero_rows_H) > 0) {
      H[non_zero_rows_H, ] <- sample_score[non_zero_rows_H, , drop = FALSE]
      has_fixed_H <- TRUE
    }
  }

  ## Iteratively estimate the NMF with Euclidean distance
  eps = 1e-12
  error_out = c()
  
  for (ite in 1:max_ite) {
    
    ## Update sample scores
    H = H * (t(W) %*% X) /
      pmax(t(W) %*% W %*% H, eps)
    
    ## Restore fixed sample scores
    if (has_fixed_H) {
      H[non_zero_rows_H, ] =
        sample_score[non_zero_rows_H, , drop = FALSE]
    }
    
    ## Enforce group structure
    component_offset = 0
    
    for (i in 1:K) {
      
      component_indices =
        (component_offset + 1):(component_offset + comp_num[i])
      
      sample_offset = 0
      
      for (j in 1:N) {
        
        sample_indices =
          (sample_offset + 1):(sample_offset + N_dataset[j])
        
        if (!(j %in% group[[i]])) {
          H[component_indices, sample_indices] = 0
        }
        
        sample_offset = sample_offset + N_dataset[j]
      }
      
      component_offset = component_offset + comp_num[i]
    }
    
    ## Update gene embeddings
    W = W * (X %*% t(H)) /
      pmax(W %*% H %*% t(H), eps)
    
    ## Restore fixed gene embeddings
    if (has_fixed_W) {
      W[, non_zero_cols_W] =
        gene_score[, non_zero_cols_W, drop = FALSE]
    }
    
    ## Calculate reconstruction error
    error_out = c(error_out, sum((X - W %*% H)^2))
    
    ## Check convergence
    if (length(error_out) >= 2) {
      
      prev_error = error_out[length(error_out) - 1]
      curr_error = error_out[length(error_out)]
      
      relative_change =
        abs(curr_error - prev_error) / max(abs(prev_error), eps)
      
      if (relative_change <= max_err) {
        break
      }
    }
  }
  

  
  
  ## Output component and scores
  list_component = list()
  list_score = list()
  for(j in 1 : N){
    list_score[[j]] = list()
  }
  
  for(i in 1 : K){
    list_component[[i]] = W[, ifelse(i == 1, 1, cumsum(comp_num)[i - 1] + 1) : cumsum(comp_num)[i], drop = FALSE]
    for(j in 1 : N){
      list_score[[j]][[i]] = H[ifelse(i == 1, 1, cumsum(comp_num)[i - 1] + 1) : cumsum(comp_num)[i], ifelse(j == 1, 1, cumsum(N_dataset)[j - 1] + 1) : cumsum(N_dataset)[j], drop = FALSE]
    }
  }
  
  
  
  ## Assign name for components
  list_component = compNameAssign(list_component, group_name)
  list_component = geneNameAssign(list_component, gene_name)
  list_score = scoreNameAssign(list_score, dataset_name, group_name)
  list_score = sampleNameAssign(list_score, sample_name)
  # list_score = filterNAValue(list_score, dataset, group)
  list_score = rebalanceData(list_score, group, dataset)
  
  ## Max-normalize sample scores while preserving reconstruction
  ## Skip components with fixed gene embeddings or sample scores
  
  eps = 1e-12
  
  fixed_components = union(non_zero_cols_W, non_zero_rows_H)
  
  component_offset = 0
  
  for (i in 1:K) {
    
    for (k in 1:comp_num[i]) {
      
      # Global component index in W and H
      global_k = component_offset + k
      
      # Skip components with fixed W or H
      if (global_k %in% fixed_components) {
        next
      }
      
      # Find maximum across all datasets participating in this group
      component_max = 0
      
      for (j in group[[i]]) {
        component_max = max(
          component_max,
          max(list_score[[j]][[i]][k, ], na.rm = TRUE)
        )
      }
      
      if (component_max > eps) {
        
        # Normalize scores across all participating datasets
        for (j in group[[i]]) {
          list_score[[j]][[i]][k, ] =
            list_score[[j]][[i]][k, ] / component_max
        }
        
        # Compensate by scaling gene embeddings
        list_component[[i]][, k] =
          list_component[[i]][, k] * component_max
      }
    }
    
    component_offset = component_offset + comp_num[i]
  }
  


  return(list(linked_component_list = list_component, score_list = list_score, error_out=error_out))
  
}
  
  