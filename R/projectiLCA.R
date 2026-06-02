#' Function to estimate sample embeddings for one dataset from a gene loading matrix derived from an iLCA analysis of another dataset.
#' 
#' projectiLCA estimates the embeddings for samples in a new dataset when given a gene loading matrix from an iLCA analysis result of another single matrix, 
#' or set of matrices (e.g. the "list_component" from a twoStageiLCA output object)
#'
#' PLEASE MAKE SURE YOUR proj_dataset AND list_component ELEMENTS HAVE MEANINGFUL ROW(GENE) NAMES - they are matched across matrices for the projection. 
#' 
#' @param proj_dataset The dataset(s) to be projected on. 
#' @param proj_group A logical vector indicating which groupings, i. e. which elements of list_component should be used for each projected dataset. The length of proj_group should match the length of list_component.
#' @param list_component a single matrix of gene loadings as a list element, or a list_component produced from a twoStageiLCA() decomposition.
#' @param ica_score ice_score produced from a twoStageiLCA() decomposition.
#' @param max_ite The maximum number of iterations for the twoStageiLCA algorithms to run, default value is set to 1000
#' @param max_err The maximum error of loss between two iterations, or the program will terminate and return, default value is set to be 0.0001
#' @param enable_normalization An argument to decide whether to use normalizaiton or not,  default is TRUE
#' @param column_sum_normalization An argument to decide whether to use column sum normalization or not, default it FALSE
#' 
#' @return A list that contains the projected scores of each dataset on every component.
#'
#' @keywords projection, twoStageiLCA, iLCA
#'
#' @examples
#' proj_dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' proj_group = c(TRUE, TRUE) # which groupings in the twoStageiLCA analysis you want to project on.
#' list_component = twoStageiLCA_res$linked_component_list # from twoStageiLCA result
#' ica_score = twoStageiLCA_res$ica_score # from twoStageiLCA result
#' res_projiLCA = projectiLCA(
#' proj_dataset = proj_dataset,
#' proj_group = proj_group,
#' list_component = list_component,
#' ica_score = ica_score)
#' 
#' 
#' 
#' @export


## Projection function

projectiLCA <- function(proj_dataset, proj_group, list_component, ica_score, max_ite = 1000, max_err = 0.0001, enable_normalization = TRUE, column_sum_normalization = FALSE){
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
    11
    
    ## Project score

    proj_sample_name = sampleNameExtractor(proj_dataset)
    group_name = names(list_component) 
    if (is.null(group_name)) {
      group_name = paste0("group", 1:length(proj_group))
    }
    proj_dataset = normalizeData(list(proj_dataset), enable_normalization, column_sum_normalization, nonnegative_normalization = TRUE)[[1]]
    

    proj_list_score = list()
    for(j in 1 : length(proj_group)){
      if(proj_group[[j]]){
        proj_list_score[[j]] = t(t((t(list_component[[j]]) %*% proj_dataset)) %*% ica_score[[j]]$K %*% ica_score[[j]]$W)
      }else{
        proj_list_score[[j]] = NA
      }
    }
    
    proj_list_score = scoreNameAssignProj(proj_list_score, group_name)
    proj_list_score = sampleNameAssignProj(proj_list_score, proj_sample_name)

    return(proj_score_list=proj_list_score)
  } else{
    stop("proj_dataset cannot be NULL")
  }
}