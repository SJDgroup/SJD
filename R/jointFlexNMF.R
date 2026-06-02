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
  error_out = c()

  for(ite in 1 : max_ite){
    H = H * (t(W) %*% X) / (t(W) %*% W %*% H)
    W = W * (X %*% t(H)) / (W %*% H %*% t(H))

    ## Restore pre-defined fixed values after each multiplicative update
    if (has_fixed_H) {
      H[non_zero_rows_H, ] <- sample_score[non_zero_rows_H, , drop = FALSE]
      ## Re-enforce group-structure zero blocks in case fixed rows overlap with them
      for(i in 1 : K){
        for(j in 1 : N){
          if(!(j %in% group[[i]])){
            H[ifelse(i == 1, 1, cumsum(comp_num)[i - 1] + 1) : cumsum(comp_num)[i],
              ifelse(j == 1, 1, cumsum(N_dataset)[j - 1] + 1) : cumsum(N_dataset)[j]] = 0
          }
        }
      }
    }
    if (has_fixed_W) {
      W[, non_zero_cols_W] <- gene_score[, non_zero_cols_W, drop = FALSE]
    }

    H[which(is.na(H))] = 0
    W[which(is.na(W))] = 0

    error_out = c(error_out, sum((X - W %*% H)^2))

    ## Break when the error difference is small
    if(length(error_out) >= 2 && abs(error_out[length(error_out)] - error_out[length(error_out) - 1]) / abs(error_out[length(error_out) - 1]) <= max_err){
      break
    }
    # print(ite)
    # print(abs(error_out[length(error_out)] - error_out[length(error_out) - 1]) / abs(error_out[length(error_out) - 1]))
  }
  
  
  # H <- t(scale(t(H), center = FALSE))
  # W = W * (X %*% t(H)) / (W %*% H %*% t(H))
  
  
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
  
  
  for(i in 1 : K){
    for(j in 1 : N){
      # Extract the subset of H based on current i and j
      list_score[[j]][[i]] <- t(apply(t(list_score[[j]][[i]]), 2, min_max_normalization))
      H[ifelse(i == 1, 1, cumsum(comp_num)[i - 1] + 1) : cumsum(comp_num)[i],
        ifelse(j == 1, 1, cumsum(N_dataset)[j - 1] + 1) : cumsum(N_dataset)[j]] <- list_score[[j]][[i]]
    }
  }
  
  #recalculate gene score
  W = W * (X %*% t(H)) / (W %*% H %*% t(H))
  
  W[is.na(W)] = 0
  
  for(i in 1 : K){
    list_component[[i]] = W[, ifelse(i == 1, 1, cumsum(comp_num)[i - 1] + 1) : cumsum(comp_num)[i], drop = FALSE]
  }

  return(list(linked_component_list = list_component, score_list = list_score, error_out=error_out))
  
}
  
  