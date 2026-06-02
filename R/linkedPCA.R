#' Linked Component Analysis
#'
#' Functions to implement the linked component analysis
#'
#' @param dataset A list of dataset to be analyzed
#' @param cov_list A list of covariance of the datasets
#' @param eigen_space A matrix of the space of the signal
#' @param group A list of grouping of the datasets, indicating the relationship between datasets
#' @param comp_num A vector indicates the dimension of each compoent
#'
#' @return A list of component
#'
#' @keywords Linked Component Analysis, sequential
#'
#' @export

linkedPCA <- function(dataset, cov_list, eigen_space, group, comp_num){
    ## Add covariance matrix together
    add.cov <- function(cov_list, group){
        add_cov_list = list()
        for(i in 1 : length(group)){
            add_cov_list[[i]] = 0
            for(j in 1 : length(cov_list)){
                if(!(j %in% group[[i]])){
                    add_cov_list[[i]] = add_cov_list[[i]] + cov_list[[j]]
                }
            }
        }
        return(add_cov_list)
    }
    ## Project the covariance onto the eigenspace
    proj.cov <- function(cov_list, eigen_space, group){
        return(lapply(cov_list, function(x){ t(eigen_space) %*% x %*% eigen_space }))
    }
    ## Find all components
    comp.all <- function(add_cov_list, comp_num, used_num){
        ## Find the component
        comp.find <- function(add_cov_list, remaining_num, comp_num, used_num){
            eigen_list = lapply(add_cov_list, eigen)
            var_list = c()
            for(i in 1 : length(add_cov_list)){
                var_list = c(var_list, ifelse(comp_num[[i]] == remaining_num[i], Inf, eigen_list[[i]]$values[length(eigen_list[[i]]$values) - used_num]))
            }

            index = which.min(var_list)
            remaining_num[index] = remaining_num[index] + 1
            temp_comp = eigen_list[[index]]$vectors[, length(eigen_list[[index]]$values) - used_num]
            used_num = used_num + 1

            add_cov_list =  lapply(add_cov_list ,function(x, temp_comp){
                x - temp_comp %*% (t(temp_comp) %*% x) - (x %*% temp_comp) %*% t(temp_comp) + temp_comp %*% (t(temp_comp) %*% x %*% temp_comp) %*% t(temp_comp)
            }, temp_comp = temp_comp)
            return(list(index = index, comp = temp_comp, remaining_num = remaining_num, add_cov_list = add_cov_list, used_num = used_num))
        }

        ## The main component of comp.all function
        remaining_num = rep(0, length(add_cov_list))
        out_comp = c()
        out_index = c()
        for(i in 1 : sum(unlist(comp_num))){
            est = comp.find(add_cov_list, remaining_num, comp_num, used_num)
            remaining_num = est$remaining_num
            add_cov_list = est$add_cov_list
            out_index = c(out_index, est$index)
            out_comp = cbind(out_comp, est$comp)
            used_num = est$used_num
        }
        return(list(out_index = out_index, out_comp = out_comp, used_num = used_num))
    }
    ## Reorganize the components into the format that can be used
    reorganize.comp <- function(cov_list, out_all, eigen_space, comp_num, group){
        out_comp <- list()
        comp_used <- rep(0, length(cov_list))
        for(i in 1 : length(group)){
            for(j in 1 : length(group[[i]])){
                comp_used[group[[i]][j]] <- comp_used[group[[i]][j]] + comp_num[[i]]
            }
            out_comp[[i]] <- vector()
        }
        for(i in 1 : length(out_all$out_index)){
            out_comp[[out_all$out_index[i]]] <- cbind(out_comp[[out_all$out_index[i]]], eigen_space %*% out_all$out_comp[, i])
        }
        for(i in 1 : length(group)){
            cov_temp = 0
            for(j in group[[i]]){
                cov_temp = cov_temp + cov_list[[j]]
            }
            out_comp[[i]] = out_comp[[i]] %*% eigen(t(out_comp[[i]]) %*% cov_temp %*% out_comp[[i]])$vectors
        }
        out_temp <- eigen_space %*% out_all$out_comp
        eigen_space <- eigen_space - out_temp %*% (t(out_temp) %*% eigen_space)
        return(list(eigen_space = eigen_space, out_comp = out_comp, comp_num = comp_num, comp_used = comp_used, used_num = out_all$used_num))
    }
    ## Output the final lca result
    lca.final <- function(dataset, cov_list, eigen_space, group, comp_num){



        ## Helper function to do the LCA decomposition
        lca.helper <- function(dataset, cov_list, eigen_space, group, comp_num, used_num){
            proj_cov_list <- proj.cov(cov_list, eigen_space, group)
            add_cov_list <- add.cov(proj_cov_list, group)
            out_all <- comp.all(add_cov_list, comp_num, used_num)
            out_reorg <- reorganize.comp(cov_list, out_all, eigen_space, comp_num, group)
            return(out_reorg)
        }

        ## Sequential algorithm to extract the components
        freq.helper = unlist(lapply(group, length))
        used_eigen_space = eigen_space
        used_num = 0
        out_comp = list()
        for(i in sort(unique(freq.helper))){
            index.helper = which(freq.helper == i)
            if(length(group[index.helper]) == 1 && length(group[index.helper][[1]]) == length(dataset)){
                out_comp[[index.helper]] = svd(used_eigen_space)$u[, 1 : comp_num[index.helper]]
            }else{
                used_lca = lca.helper(dataset, cov_list, used_eigen_space, group[index.helper], comp_num[index.helper], used_num)
                used_eigen_space = used_lca$eigen_space
                used_num = used_lca$used_num
                out_comp[index.helper] = used_lca$out_comp
            }
        }

        ## Parameters to be initialized
        N = length(dataset)
        K = length(group)
        M = sum(comp_num)
        p = nrow(dataset[[1]])
        N_dataset = unlist(lapply(dataset, ncol))

        ## Output the component and scores
        list_score = list()
        list_component = list()

        for(j in 1 : N){
            list_score[[j]] = list()
        }

        for(i in 1 : K){
            list_component[[i]] = matrix(0, nrow = p, ncol = comp_num[i])
            for(j in 1 : N){
                list_score[[j]][[i]] = matrix(0, nrow = comp_num[i], ncol = N_dataset[j])
            }
        }


        ## output the list component
        for(i in 1 : K){
            list_component[[i]] = out_comp[[i]]
            for(j in 1 : N){
                if(j %in% group[[i]]){
                    list_score[[j]][[i]] = t(list_component[[i]]) %*% dataset[[j]]
                }
            }
        }


        return(list(linked_component_list = list_component, score_list = list_score))
    }

    return(lca.final(dataset, cov_list, eigen_space, group, comp_num))
}
