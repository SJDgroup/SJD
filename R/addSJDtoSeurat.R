#' Add SJD loadings to a Seurat object
#'
#' @param Seurat.obj A Seurat Object.
#' @param SJDoutput A SJD object of decomposition results.
#' @param Dataset The name of the Dataset of interests in the SJD score_list component
#' @param SJDloading The name of the SJDloading of interests in the SJD score_list component under Dataset.
#' @param SJDmethod The name of the SJD method used for the decomposition, e.g. "twoStageLCA". If NA, will use "SJD"
#'
#' @export

addSJDtoSeurat = function(Seurat.obj,SJDoutput,Dataset,SJDloading,SJDmethod=NA){

    if(is.na(SJDmethod)){key = 'SJD_'}
    else{key = paste0(SJDmethod,'_')}

    embeddings = t(SJDoutput$score_list[[Dataset]][[SJDloading]])
    rownames(embeddings)=paste0(SJDloading,'_PVE',gsub('.*: ','',rownames(embeddings)))
    loadings = SJDoutput$linked_component_list[[SJDloading]]

    stopifnot("SJD dimension must match Seurat dimension!" =
                  identical(colnames(embeddings),colnames(seurat.obj)) &
                  NCOL(embeddings)==NCOL(seurat.obj))

    Seurat.obj = CreateDimReducObject(embeddings = embeddings, loadings = loadings, key = key,
                                      assay = DefaultAssay(seurat.obj))

    return(Seurat.obj)
}
