#' Dataset Name Extractor
#'
#' Extract data.frame name from a list of data.frames, if
#' no name is given in the list, it will be renamed as
#' "dataset_no.index"
#'
#' @param dataset A list of data sets
#'
#' @return A vector of strings of dataset name
#'
#' @keywords name
#'
#' @examples
#' dataset = list(
#' x = matrix(c(1 : 4), nrow = 2),
#' y = matrix(c(1 : 4), nrow = 2))
#'
#' datasetNameExtractor(dataset)
#'
#' @export

datasetNameExtractor <- function(dataset){
    out = names(dataset)
    if(is.list(dataset)){
      if(is.null(out)){
        out = rep("", length(dataset))
      }
    }else{
      out = ""
    }
    index = 1
    for(i in 1 : length(out)){
        if(out[i] == ""){
            out[i] = paste0("dataset_No.", i)
        }
    }
    return(out)
}



#' Group Name Extractor
#'
#' Extract group name from groups, if no name is given in the input
#' list, the output vector will be represented as
#' "component_No.index"
#'
#' @param group A list of group assignments for the data sets with group name on it
#'
#' @return A vector of strings of group name
#'
#' @keywords name
#'
#' @examples
#' dataset = list(
#' x = c(1, 2, 3),
#' y = c(1, 2, 4))
#'
#' @export

groupNameExtractor <- function(group){
    out = names(group)
    if(is.null(out)){
        out = rep("", length(group))
    }
    index = 1
    for(i in 1 : length(out)){
        if(out[i] == ""){
            out[i] = paste0("component_No.", i)
        }
    }
    return(out)
}



#' Gene Name Extractor
#'
#' Extract gene names from the input dataset
#'
#' @param dataset A list of data sets to be analyzed
#'
#' @return a vector of strings of gene name
#'
#' @keywords gene, name
#'
#' @examples
#' x = matrix(c(1 : 4), nrow = 2)
#' rownames(x) = c("row1", "row2")
#' y = matrix(c(1 : 4), nrow = 2)
#' rownames(y) = c("row1", "row2")
#' dataset = list(x, y)
#' geneNameExtractor(dataset)
#'
#' @export

geneNameExtractor <- function(dataset){
    gene_name = rownames(dataset[[1]])
    if(is.null(gene_name)){
        gene_name = rep("", nrow(dataset[[1]]))
    }

    for(i in 1 : nrow(dataset[[1]])){
        if(gene_name[i] == ""){
            gene_name[i] = paste0("gene_No.", i)
        }
    }

    return(gene_name)
}


#' Sample name extractor
#'
#' Extract a list of sample names from input list of datasets
#'
#' @param dataset A list of datasets containing sample names
#'
#' @return A vector of sample names
#'
#' @keywords name, extractor
#'
#' @examples
#' x = matrix(c(1 : 4), nrow = 2)
#' colnames(x) = c("sp1", "sp2")
#' y = matrix(c(1 : 4), nrow = 2)
#' colnames(y) = c("sp3", "sp4")
#' dataset = list(x, y)
#' sampleNameExtractor(dataset)
#'
#' @export

sampleNameExtractor <- function(dataset){
    sample_name = list()
    if(is.list(dataset)){
      for(i in 1 : length(dataset)){
        temp_name = colnames(dataset[[i]])
        if(is.null(temp_name)){
          sample_name[[i]] = rep("", ncol(dataset[[i]]))
        }else{
          sample_name[[i]] = temp_name
        }
        for(j in 1 : length(sample_name[[i]])){
          if(sample_name[[i]][j] == ""){
            sample_name[[i]][j] = paste0("subject_No.", j)
          }
        }
      }
    }else{
      temp_name = colnames(dataset)
      if(is.null(temp_name)){
        sample_name = rep("", ncol(dataset))
      }else{
        sample_name = temp_name
      }
      for(j in 1 : length(sample_name)){
        if(sample_name[j] == ""){
          sample_name[j] = paste0("subject_No.", j)
        }
      }
    }

    return(sample_name)
}
