#' genetrate new UMAP and TSNE coordinates given LCA object
#'
#' @param twoStageLCAobj list output of tsLCA run, list
#' @param group name of the weights group from the SJD_algorithm output, i.e 'Shared.All.13'
#' @param n_comp number of components to make umaps and tsne, int
#' @param add_to_meta adds tsne and umap coordinates to metadata, bool
#' @param meta_list Meta.List that was passed into SJD run, str
#' @export

umap_tsne_onLC=function(twoStageLCAobj,group,n_comp,add_to_meta,meta_list){


  ## add argument to add umap/tsne coords to metadata obj
  for (i in names(twoStageLCAobj$score_list)){
    scores_sample_shared=t(twoStageLCAobj$score_list[[i]][[group]])
    scores_sample_shared=as.data.frame(scores_sample_shared)
    shared_umap=uwot::umap(scores_sample_shared[,1:n_comp], n_neighbors = 50, learning_rate = 0.5, init = "random")
    shared_tsne=Rtsne::Rtsne(scores_sample_shared[,1:n_comp])
    umaptsne_df=as.data.frame(shared_umap[,1])
    colnames(umaptsne_df)=c('shared_UMAP_1')
    umaptsne_df['shared_UMAP_2']=shared_umap[,2]
    umaptsne_df['shared_tSNE_1']=shared_tsne$Y[,1]
    umaptsne_df['shared_tSNE_2']=shared_tsne$Y[,2]

    if (add_to_meta==TRUE){
      method=paste(n_comp,'LCs',sep='')
      meta_list[[i]][paste('tsLCA',method,'UMAP_1',sep = '_')]=shared_umap[,1]
      meta_list[[i]][paste('tsLCA',method,'UMAP_2',sep = '_')]=shared_umap[,2]
      meta_list[[i]][paste('tsLCA',method,'tSNE_1',sep = '_')]=shared_tsne$Y[,1]
      meta_list[[i]][paste('tsLCA',method,'tSNE_2',sep = '_')]=shared_tsne$Y[,2]
    }

  }
if (add_to_meta==TRUE){results=list(shared_umap1=shared_umap[,1],shared_umap2=shared_umap[,2],shared_tsne1=shared_tsne$Y[,1],shared_tsne2=shared_tsne$Y[,2],meta_list=meta_list)}#note from carlo: use this if dropping plotting
if (add_to_meta==FALSE){results=list(shared_umap1=shared_umap[,1],shared_umap2=shared_umap[,2],shared_tsne1=shared_tsne$Y[,1],shared_tsne2=shared_tsne$Y[,2])}#note from carlo: use this if dropping plotting
  return (results)
}
