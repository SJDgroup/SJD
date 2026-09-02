#' Function to estimate sample embeddings for one dataset from a gene loading matrix derived from an NMF decomposition of another dataset.
#' 
#' projectNMF estimates the embeddings for samples in a new dataset when given a gene loading matrix from an NMF decomposition of another single matrix, 
#' or set of matrices (e.g. the "list_component" from a jointNMF output object)
#'
#' PLEASE MAKE SURE YOUR proj_dataset AND list_component ELEMENTS HAVE MEANINGFUL ROW(GENE) NAMES - they are matched across matrices for the projection. 
#' 
#' list_component should be jointNMF$linked_component_list from the jointNMF result
#' 
#' @param proj_dataset The dataset to be projected on. 
#' @param proj_group A logical vector indicating which groupings, i. e. which elements of list_component should be used for each projected dataset. The length of proj_group should match the length of list_component.
#' @param list_component a single matrix of gene loadings as a list element, or a list_component produced from a jointNMF() decomposition.
#' @param max_ite The maximum number of iterations for the jointNMF algorithms to run, default value is set to 1000
#' @param max_err The maximum error of loss between two iterations, or the program will terminate and return, default value is set to be 0.0001
#' @param enable_normalization An argument to decide whether to use normalizaiton or not,  default is TRUE
#' @param column_sum_normalization An argument to decide whether to use column sum normalization or not, default it FALSE
#' 
#' @return A list that contains the 1] projected scores of each dataset on every component. and 2] the log of errors as the NMF was iterated.
#'
#' @keywords projection, joint, NMF
#'
#' @examples
#' proj_dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' proj_group = c(TRUE, TRUE)
#' list_component = list(matrix(runif(600, 1, 2), nrow = 100, ncol = 6), matrix(runif(200, 1, 2), nrow = 100, ncol = 2))
#' res_projNMF = projectNMF(
#' proj_dataset = proj_dataset,
#' proj_group = proj_group,
#' list_component = list_component)
#' 
#' 
#'
#' @export


## Projection function

projectNMF <- function(proj_dataset, proj_group, list_component, max_ite = 1000, max_err = 0.0001, enable_normalization = TRUE, column_sum_normalization = FALSE){
  if(!is.null(proj_dataset)){
    if (length(proj_group) != length(list_component)){
      stop("Error:length of proj_group should equal to length of list_component.")
    }
    
    proj_dataset = as.matrix(proj_dataset)
      
    if(!is.null(rownames(proj_dataset))){
      CMNgenes=rownames(proj_dataset)[rownames(proj_dataset)%in%rownames(list_component[[1]])]
      common_genes.DATA = match(CMNgenes,rownames(proj_dataset))
      common_genes.COMP = match(CMNgenes,rownames(list_component[[1]]))
      orig_gene_length = nrow(proj_dataset)
      proj_dataset = proj_dataset[common_genes.DATA,]
      for (i in 1:length(list_component)) {
        list_component[[i]] = list_component[[i]][common_genes.COMP,,drop=FALSE]
      }
      print(paste("Input", orig_gene_length, "genes in proj_dataset, found", nrow(proj_dataset), "genes in common."))
    }
    
    
    p = nrow(proj_dataset)
    comp_num = unlist(lapply(list_component, function(x) if (length(x) > 0) ncol(x)))
    
    # Initialize the matrix W with zeros
    W <- matrix(0, nrow = p, ncol = sum(comp_num))
    K = length(proj_group)
    # Reconstruct matrix W using list_component and comp_num
    start_col <- 1
    for (i in 1:K) {
      end_col <- start_col + comp_num[i] - 1
      W[, start_col:end_col] <- list_component[[i]]
      start_col <- end_col + 1
    }
  
    
    proj_sample_name = sampleNameExtractor(proj_dataset)
    group_name = names(list_component) 
    if (is.null(group_name)) {
      group_name = paste0("group", 1:length(proj_group))
    }
    proj_dataset = normalizeData(list(proj_dataset), enable_normalization, column_sum_normalization, nonnegative_normalization = TRUE)[[1]]
    
    M = sum(comp_num) 
    p = nrow(proj_dataset)
    col = ncol(proj_dataset)
    
    
    ## Initialize the W and H for Nonnegative Matrix Factorization
    max_element = -Inf
    min_element = Inf
    max_element = max(max_element, max(proj_dataset))
    min_element = min(min_element, min(proj_dataset))
    
    
    ## Initialize the values of W and H
    X = proj_dataset
    
    
    H = c()
    for(i in 1 : K){ # K= length(group) = 3 in the example
      H_temp = c()
      H_temp = matrix(runif(col * comp_num[i], min_element, max_element), nrow = comp_num[i], ncol = col)
      H = rbind(H, H_temp)
    }
    
    
    ## Iteratively estimate the NMF with Euclidean distance
    error_out = c()
    
    for(ite in 1 : max_ite){
      ## H is score (sample matrix S)
      H = H * (t(W) %*% X) / (t(W) %*% W %*% H)# W is list_component (common gene matrix G)
      
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
        if(comp_num[[j]] == 1) {
          proj_list_score[[j]] = matrix(proj_list_score[[j]], nrow=1, ncol=col)
        }
      } else {
        proj_list_score[[j]] = matrix(0, nrow = comp_num[j], ncol = col)
      }
    }

    proj_list_score = scoreNameAssignProj(proj_list_score, group_name)
    proj_list_score = sampleNameAssignProj(proj_list_score, proj_sample_name)
    for(i in 1 : length(proj_group)){
      if(proj_group[i]){
        proj_list_score[[i]] = proj_list_score[[i]] * sqrt(ncol(proj_dataset))
      }
    }
    return(list(proj_score_list=proj_list_score,error_out=error_out))
  } else{
    stop("proj_dataset cannot be NULL")
  }
}





