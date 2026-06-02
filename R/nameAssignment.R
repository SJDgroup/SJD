#' Score Name Assign for Seperate Analysis
#'
#' Assign name to scores based on the extracted name for each dataset,
#' it takes the vector of dataset_name outputed by  the 'datasetNameExtractor'
#' function
#'
#' @param score_list List of scores in sep analysis
#' @param dataset_name List of dataset names extracted by datasetNameExtractor
#'
#' @return A list of scores, content same as input name changed based on dataset_name
#'
#' @keywords score, name
#'
#' @examples
#' score_list = list(matrix(c(1 : 4), nrow = 2), matrix(c(5 : 8), nrow = 2))
#' dataset_name = c("x", "y")
#' scoreNameAssignSep(score_list, dataset_name)
#'
#' @export
#' 
scoreNameAssignSep <- function(score_list, dataset_name){
    for(i in 1 : length(score_list)){
        names(score_list)[i] = dataset_name[i]
        rownames(score_list[[i]]) = sapply(1 : nrow(score_list[[i]]), FUN = function(x){paste0(dataset_name[i], "_", "subcomp.", x)} )
    }
    return(score_list)
}

#' Score Name Assignment for projectNMF
#'
#' Assign names to score lists based on group_name,
#' the group_name
#' is the output of the groupNameExtractor
#'
#' @param score_list A list of scores in the analysis
#' @param group_name A vector of group names extracted by groupNameExtractor function
#'
#' @return A list of scores assigned with names
#'
#' @keywords score, name
#'
#' @examples
#' score_list = list(
#' list(matrix(c(1 : 4), nrow = 2), matrix(c(1 : 4), nrow = 2)),
#' list(matrix(c(1 : 4), nrow = 2), matrix(c(1 : 4), nrow = 2)))
#' group_name = c("comp1", "comp2")
#' scoreNameAssignProj(score_list, group_name)
#'
#' @export

scoreNameAssignProj <- function(score_list, group_name){
  for(i in 1 : length(score_list)){
    names(score_list)[i] = group_name[i]
    rownames(score_list[[i]]) = sapply(1 : nrow(score_list[[i]]), FUN = function(x){paste0(group_name[i], "_", "subcomp.", x)} )
  }
  return(score_list)
}

#' Score Name Assignment for Concatenate, Joint and TwoStageLCA Analysis
#'
#' Assign names to score lists based on dataset_name and group_name,
#' the dataset_name is the output of the dataNameExtractor, the group_name
#' is the output of the groupNameExtractor
#'
#' @param score_list A list of scores in the analysis
#' @param dataset_name A vector of datasets names extracted by dataNameExtractor function
#' @param group_name A vector of group names extracted by groupNameExtractor function
#'
#' @return A list of scores assigned with names
#'
#' @keywords score, name
#'
#' @examples
#' score_list = list(
#' list(matrix(c(1 : 4), nrow = 2), matrix(c(1 : 4), nrow = 2)),
#' list(matrix(c(1 : 4), nrow = 2), matrix(c(1 : 4), nrow = 2)))
#' dataset_name = c("dat1", "dat2")
#' group_name = c("comp1", "comp2")
#' scoreNameAssign(score_list, dataset_name, group_name)
#'
#' @export

scoreNameAssign <- function(score_list, dataset_name, group_name){
    for(i in 1 : length(score_list)){
        names(score_list)[i] = dataset_name[i]
        for(j in 1 : length(score_list[[i]])){
            names(score_list[[i]])[j] = group_name[j]
            if(!any(is.na(score_list[[i]][[j]]))){
                rownames(score_list[[i]][[j]]) = sapply(1 : nrow(score_list[[i]][[j]]), FUN = function(x){paste0(group_name[j], "_", "subcomp.", x)} )
            }
        }
    }
    return(score_list)
}

#' Components Name Assignment for Seperate Analysis
#'
#' Assign name to components in the seperate analysis, sepPCA, sepICA, sepNMF
#'
#' @param linked_component_list list of components extracted
#' @param dataset_name A vector of names for the datasets
#'
#' @return renamed list of linked_component_list
#'
#' @keywords component, name
#'
#' @examples
#' linked_component_list = list(matrix(c(1:4), nrow = 2), matrix(c(1:4), nrow = 2))
#' dataset_name = c("x", "y")
#' compNameAssignSep(linked_component_list, dataset_name)
#'
#' @export

compNameAssignSep <- function(linked_component_list, dataset_name){
    for(i in 1 : length(linked_component_list)){
        names(linked_component_list)[i] = dataset_name[i]
        colnames(linked_component_list[[i]]) = sapply(1 : ncol(linked_component_list[[i]]), FUN = function(x){paste0(dataset_name[i], "_", "subcomp.", x)} )
    }
    return(linked_component_list)

}

#' Components Name Assignment for concatenate, joint and twoStageLCA analysis
#'
#' Assign name to components in the concatenate, joint and twoStageLCA methods
#'
#' @param linked_component_list list of components extracted
#' @param group_name A vector of names for the datasets
#'
#' @return renamed list of linked_component_list
#'
#' @keywords component, name
#'
#' @examples
#' linked_component_list = list(matrix(c(1:4), nrow = 2), matrix(c(1:4), nrow = 2))
#' group_name = c("x", "y")
#' compNameAssign(linked_component_list, group_name)
#'
#' @export

compNameAssign <- function(linked_component_list, group_name){
    for(i in 1 : length(linked_component_list)){
        names(linked_component_list)[i] = group_name[i]
        colnames(linked_component_list[[i]]) = sapply(1 : ncol(linked_component_list[[i]]), FUN = function(x){paste0(group_name[i], "_", "subcomp.", x)} )
    }
    return(linked_component_list)
}

#' Gene name assignment
#'
#' Assign gene names to component derived in the analysis
#'
#' @param linked_component_list A list of extracted components
#' @param gene_name A vector of strings for gene names
#'
#' @return A list of components with gene name added
#'
#' @keywords gene, name
#'
#' @examples
#' x = matrix(c(1 : 4), nrow = 2)
#' y = matrix(c(1 : 4), nrow = 2)
#' component_list = list(x, y)
#' gene_name = c("gene.1", "gene.2")
#' geneNameAssign(component_list, gene_name)
#'
#' @export

geneNameAssign <- function(linked_component_list, gene_name){
    for(i in 1 : length(linked_component_list)){
        rownames(linked_component_list[[i]]) = gene_name
    }
    return(linked_component_list)
}



#' Sample Name Assignment for projectNMF
#'
#' Assign sample names to dataset,
#' it takes two arguments, the first is the list of scores,
#' the second is the list of sample names
#'
#' @param score_list List of score for each data.frame
#' @param sample_name List of names for the samples in the list
#'
#' @return A list of scores for the samples
#'
#' @keywords sample, name
#'
#' @examples
#' x = matrix(c(1:4), nrow = 2)
#' y = matrix(c(1:4), nrow = 2)
#' score_list = list(x, y)
#' sample_name = list("x.sample.1", "x.sample.2")
#' sampleNameAssignProj(score_list, sample_name)
#'
#' @export


sampleNameAssignProj <- function(score_list, sample_name){
  for(i in 1 : length(score_list)){
    colnames(score_list[[i]]) = sample_name
  }
  return(score_list)
}

#' Sample Name Assignment for Seperate Analysis
#'
#' Assign sample names to each dataset in the list of data sets,
#' it takes two arguments, the first is the list of scores,
#' the second is the list of sample names
#'
#' @param score_list List of score for each data.frame
#' @param sample_name List of names for the samples in the list
#'
#' @return A list of scores for the samples
#'
#' @keywords sample, name
#'
#' @examples
#' x = matrix(c(1:4), nrow = 2)
#' y = matrix(c(1:4), nrow = 2)
#' score_list = list(x, y)
#' sample_name = list(c("x.sample.1", "x.sample.2"), c("y.sample.1", "y.sample.2"))
#' sampleNameAssignSep(score_list, sample_name)
#'
#' @export

sampleNameAssignSep <- function(score_list, sample_name){
    for(i in 1 : length(score_list)){
        colnames(score_list[[i]]) = sample_name[[i]]
    }
    return(score_list)
}

#' Sample Name Assignment for Concatenate, Joint and TwoStageLCA Analysis
#'
#' Assign sample names to the scores
#'
#' @param score_list A list of scores to do analysis
#' @param sample_name A list of names for samples to be analyzed
#'
#' @return A list of scores to analyze
#'
#' @keywords sample, name
#'
#' @examples
#' x = list(matrix(c(1 : 4), nrow = 2), matrix(c(1 : 4), nrow = 2))
#' y = list(matrix(c(1 : 4), nrow = 2), matrix(c(1 : 4), nrow = 2))
#' score_list = list(x, y)
#' sample_name = list(c("x.sample.1", "x.sample.2"), c("y.sample.1", "y.sample.2"))
#' sampleNameAssign(score_list, sample_name)
#'
#' @export

sampleNameAssign <- function(score_list, sample_name){
    for(i in 1 : length(score_list)){
        for(j in 1 : length(score_list[[i]])){
            if(!any(is.na(score_list[[i]][[j]]))){
                colnames(score_list[[i]][[j]]) = sample_name[[i]]
            }
        }
    }
    return(score_list)
}
