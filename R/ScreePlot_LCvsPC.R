#' genetrate scree plot using tsLCA and sepPCA embeddings for each matrix
#'
#' @param LCAobj list output of tsLCA run, list
#' @param group name of the weights group from the SJD_algorithm output, i.e 'Shared.All.13', str
#' @param sepPCAobj list output of sepPCA run, list
#' @export

ScreePlot_LCvsPC=function(LCAobj,group,sepPCAobj){
  dims2=1:length(colnames(LCAobj$linked_component_list[[group]]))
  pve_df=data.frame(dims2)
  dataset=c()
  j=1
  for (i in LCAobj$score_list){
    pve=c()
    for (k in rownames(i[[group]])){
      pve=c(pve,as.double(strsplit(k,' ')[[1]][3]))
    }
    pve_df[paste('tsLCA',names(LCAobj$score_list[j]),sep='_')]=pve
    j=j+1
  }
  ### get pve's from sepPCAobj object ###
  dataset=c()
  j=1
  for (i in sepPCAobj$score_list){
    pve=c()
    for (k in rownames(i)){
      pve=c(pve,as.double(strsplit(k,' ')[[1]][3]))
    }
    pve_df[paste('PCA',names(sepPCAobj$score_list[j]),sep='_')]=pve
    j=j+1
  }
  n <- length(names(LCAobj$score_list))
  return(pve_df)
}
