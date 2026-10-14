#' Flexible Decomposition with Nonnegative Matrix Factorization
#' 
#' Function to decompose a matrix with known gene loading and/or sample embedding with Nonnegative Matrix Factorization.
#' 
#' flexNMF decomposes the matrix while keeps the given/known gene columns or sample rows fixed, estimating the rest of the gene and sample matrices. 
#'
#'
#' dataset, gene_score=NULL, sample_score=NULL, normalize = TRUE, max_ite = 1000, max_err = 0.0001, enable_normalization = TRUE, column_sum_normalization = FALSE
#' @param dataset The dataset to be decomposed 
#' @param gene_score The known gene embedding columns/patterns to be fixed
#' @param sample_score The known sample score columns/patterns to be fixed
#' @param normalize whether the dataset matrix needs to be normalized
#' @param max_ite The maximum number of iterations for the jointNMF algorithms to run, default value is set to 1000
#' @param max_err The maximum error of loss between two iterations, or the program will terminate and return, default value is set to be 0.0001
#' @param enable_normalization An argument to decide whether to use normalizaiton or not,  default is TRUE
#' @param column_sum_normalization An argument to decide whether to use column sum normalization or not, default it FALSE
#' 
#' @return A list that contains the 1] gene embedding scores of the decomposed dataset. 2] sample scores of the decmposed dataset. and 3] the log of errors as the NMF was iterated.
#'
#' @keywords flexible, NMF, decomposition
#'
#' @examples
#' dataset = matrix(runif(5000, 1, 2), nrow = 100, ncol = 50)
#' gene_score = matrix(runif(200, 1, 2), nrow = 100, ncol = 2)
#' sample_score = matrix(runif(100, 1, 2), nrow = 2, ncol = 50)
#' normalize = TRUE
#' res_flexNMF = flexNMF(
#' dataset = dataset,
#' gene_score = gene_score,
#' sample_score = sample_score,
#' normalize = normalize)
#' 
#'
#' @export



flexNMF = function(dataset, gene_score=NULL, sample_score=NULL, normalize = TRUE, max_ite = 1000, max_err = 0.0001, enable_normalization = TRUE, column_sum_normalization = FALSE){
  if(!is.null(dataset)){ 
    if(!is.null(gene_score) && !is.null(sample_score)){
      
      num_patterns_g = ncol(gene_score)
      num_patterns_s = nrow(sample_score)
      
      if(num_patterns_g != num_patterns_s){
        stop("ncol(gene_score) must match nrow(sample_socre)")
      }
      
      dataset = as.matrix(dataset)
      
      p = nrow(dataset)
      col = ncol(dataset)
      
      non_zero_cols_index = which(colSums(gene_score) > 0)
      non_zero_cols = gene_score[, non_zero_cols_index,drop = FALSE]
      
      non_zero_rows_index = which(rowSums(sample_score) > 0)
      non_zero_rows = sample_score[non_zero_rows_index, , drop = FALSE]
      
      
      ## Initialize the W and H for Nonnegative Matrix Factorization
      max_element = -Inf
      min_element = Inf
      max_element = max(max_element, max(dataset))
      min_element = min(min_element, min(dataset))
      
      # Initialize the matrix W with random values (gene score)
      W = matrix(runif(p * num_patterns_g, min(non_zero_cols), max(non_zero_cols)), nrow = p, ncol = num_patterns_g)
      # Initialize the matrix H with random values (sample score)
      H = matrix(runif(num_patterns_g * col, min(non_zero_rows), max(non_zero_rows)), nrow = num_patterns_g, ncol = col)
      
      #extract non-zero columns and rows
      W[, non_zero_cols_index] = non_zero_cols
      H[non_zero_rows_index, ] = non_zero_rows
      
    
      # create masks for H and W
      mask_W = matrix(1, nrow = p, ncol = num_patterns_g)
      mask_W[, non_zero_cols_index] = 0
      mask_H = matrix(1, nrow = num_patterns_g, ncol= col)
      mask_H[non_zero_rows_index, ] = 0
      
      
      proj_sample_name = sampleNameExtractor(dataset)
      
      if (normalize){
        dataset = normalizeData(list(dataset), enable_normalization, column_sum_normalization, nonnegative_normalization = TRUE)[[1]]
      }
      
      X = dataset
      
      ## Iteratively estimate the NMF with Euclidean distance
      error_out = c()
      
      W_temp = W
      H_temp = H
      
      eps = 1e-12
      
      for (ite in 1:max_ite) {
        
        H_temp = H_temp * (t(W_temp) %*% X) /
          pmax(t(W_temp) %*% W_temp %*% H_temp, eps)
        
        H_temp[non_zero_rows_index, ] = non_zero_rows
        
        W_temp = W_temp * (X %*% t(H_temp)) /
          pmax(W_temp %*% H_temp %*% t(H_temp), eps)
        
        W_temp[, non_zero_cols_index] = non_zero_cols
        
        error_out = c(error_out, sum((X - W_temp %*% H_temp)^2))
        
        if (length(error_out) >= 2) {
          
          prev_error = error_out[length(error_out) - 1]
          curr_error = error_out[length(error_out)]
          
          relative_change = abs(curr_error - prev_error) /
            max(abs(prev_error), eps)
          
          if (relative_change <= max_err) {
            break
          }
        }
      }
      
      W = W_temp
      H = H_temp
      
      ## Max-normalize only components without fixed constraints
      eps = 1e-12
      
      fixed_components = union(
        non_zero_cols_index,
        non_zero_rows_index
      )
      
      for (k in 1:nrow(H)) {
        
        # Do not modify components with fixed W or fixed H
        if (k %in% fixed_components) {
          next
        }
        
        component_max = max(H[k, ], na.rm = TRUE)
        
        if (component_max > eps) {
          
          H[k, ] = H[k, ] / component_max
          
          W[, k] = W[, k] * component_max
        }
      }
      
      return(list(res_gene_score=W, res_sample_score=H, error_out=error_out))
      
    } else{
      stop("gene_score and sample_score cannot be NULL")
    }
  } else{
    stop("dataset cannot be NULL")
  }
}
