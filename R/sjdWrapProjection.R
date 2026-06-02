#' SJD Wrap for projection
#'
#' Wrapping up a new list of projection matrices with the list computed from the function of sjdWrap
#'
#' @param data.list.template input list of expression matrices from the output of the sjdWrap function
#' @param data.list.projection input list of expression matrices needed to be matched with data.list.template
#' @param species.template character of species type from the sjdWrap function
#' @param species.vec.projection list of species of each matrix to be projected
#' @param geneType.template character of gene/rowname type from the sjdWrap function
#' @param geneType.vec.projection list of output gene/rowname type of each matrix to be projected
#'
#' @importFrom dplyr inner_join
#'
#' @return A list of expression matrices (of different species) with only shared genes to be projected
#'
#' @keywords shared gene expression
#'
#' @examples
#'
#' ## Load NeuroGenesis4 data into R
#' data(NeuroGenesis4)
#'
#' ## sjdWrap of the training data sets
#' SJDdataIN = sjdWrap(
#' data.list = NeuroGenesis4,
#' species.vector=c("human","human","human","mouse"),
#' geneType.vector=c("symbol","ensembl","symbol","symbol"),
#' geneType.out="symbol",
#' species.out="human")
#'
#' ## Sample from data, serving as the projection expression matrices
#' NeuroGenesis4.sample = NeuroGenesis4
#' NeuroGenesis4.sample[[1]] = NeuroGenesis4.sample[[1]][-5,]
#' rownames(NeuroGenesis4.sample[[1]])[5] = paste0(rownames(NeuroGenesis4.sample[[1]])[5], ".test")
#'
#' NeuroGenesis4.sample[[2]] = NeuroGenesis4.sample[[2]][-10,]
#' rownames(NeuroGenesis4.sample[[2]])[10] = paste0(rownames(NeuroGenesis4.sample[[2]])[10], ".test")
#'
#' NeuroGenesis4.sample[[3]] = NeuroGenesis4.sample[[3]][-15,]
#' rownames(NeuroGenesis4.sample[[3]])[15] = paste0(rownames(NeuroGenesis4.sample[[3]])[15], ".test")
#'
#' NeuroGenesis4.sample[[4]] = NeuroGenesis4.sample[[4]][-20,]
#' rownames(NeuroGenesis4.sample[[4]])[20] = paste0(rownames(NeuroGenesis4.sample[[4]])[20], ".test")
#'
#' SJDdataProjection = sjdWrapProjection(
#' SJDdataIN, NeuroGenesis4.sample, "human", c("human","human","human","mouse"),
#' "symbol", c("symbol","ensembl","symbol","symbol"))
#'
#' @export

sjdWrapProjection <- function(data.list.template, data.list.projection,
                              species.template, species.vec.projection,
                              geneType.template, geneType.vec.projection
                              ){
    N = length(data.list.projection) #number of datasets

    if(N != length(species.vec.projection)){
        print("length of data.list != length(species.vector)")
        return(NULL)
    }

    if(N != length(geneType.vec.projection)){
        print("length of data.list != length(geneType.vector)")
        return(NULL)
    }


    genes.tbl.list = list()
    new.data.list = list()


    cat('Using biomaRt to connect gene IDs across', N, 'datasets:\n')

    for (i in 1 : N) {

        cat('Getting biomaRt IDs for dataset',i,'\n')

        species.1 = species.vec.projection[i]
        genes.1 = rownames(data.list.projection[[i]])
        geneType.1 = geneType.vec.projection[i]

        genes.tbl = getMatch(genes.1, inSpecies = species.1, inType = geneType.1, newSpecies = species.template)
        genes.tbl = na.omit(genes.tbl)

        genes.tbl$ORIGid = genes.tbl[, "genes"]

        if(geneType.template == "symbol") {
            genes.tbl[, "genes"] = genes.tbl[, 5]
        } # geneType.out

        if(geneType.template == "ensembl"){
            genes.tbl[, "genes"] = genes.tbl[, 6]
        } # geneType.out

        indxNoDup = !duplicated(genes.tbl[, "genes"]) & !is.na(genes.tbl[, "genes"]) & genes.tbl[, "genes"] != ""
        genes.tbl.list[[i]] = genes.tbl[indxNoDup, c("genes", "ORIGid")]
        rm(genes.tbl)
    }


    ## Extract the intersection between genes.tbl.list genes and template genes
    cat('constructed', length(genes.tbl.list), 'tables of cross-species matching genes\n')

    genes.template.vec = rownames(data.list.template[[1]])
    for(i in 1 : N){
        genes.tbl.list[[i]] = genes.tbl.list[[i]][match(intersect(genes.template.vec, genes.tbl.list[[i]]$genes), genes.tbl.list[[i]]$genes), ]
        new.data.list[[i]] = data.frame(matrix(0, nrow = nrow(data.list.template[[1]]), ncol = ncol(data.list.projection[[i]])))

        new.data.list[[i]][match(genes.tbl.list[[i]]$genes, genes.template.vec), ] = data.list.projection[[i]][match(genes.tbl.list[[i]]$ORIGid, rownames(data.list.projection[[i]])), ]

        rownames(new.data.list[[i]]) = genes.template.vec
        colnames(new.data.list[[i]]) = colnames(data.list.projection[[i]])
    }

    #
    # genes.tbl = inner_join(genes.tbl.list[[1]], genes.tbl.list[[2]], by = 'genes', na_matches = "never") #geneType.out
    #
    # if(N >= 3){
    #     for(i in 3 : N){
    #         genes.tbl = inner_join(genes.tbl, genes.tbl.list[[i]], by = 'genes', na_matches = "never") #geneType.out
    #     }
    # }
    #
    # cat('we found', dim(genes.tbl)[1], 'shared genes in', length(genes.tbl.list), 'datasets\n')
    #
    # for (i in 1 : N) {
    #     indxR = match(genes.tbl[, i + 1], rownames(data.list[[i]])) #ORIGid
    #     data = data.list[[i]][indxR, ]
    #     rownames(data) = genes.tbl[, "genes"]
    #     new.data.list[[i]] = as.data.frame(data)
    #     rm(data)
    # }
    #
    # cat('new data list of', length(new.data.list), 'datasets constructed\n')

    names(new.data.list) = names(data.list.projection)

    return(new.data.list)
}
