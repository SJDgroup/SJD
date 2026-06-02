#' SJD Wrap
#'
#' wrapping up expression matrices (of different species) with only shared genes as SJD input
#'
#' @param data.list input list of expression matrices from different species i.e human and mouse datasets
#' @param species.vector character of species type of each matrix i.e c('human', 'mouse', 'mouse')
#' @param geneType.vector character of gene/rowname type of each matrix i.e c("symbol","ensembl","symbol")
#' @param geneType.out character of output gene/rowname type of each matrix i.e "symbol"
#' @param species.out character of output species type for gene/rowname
#'
#' @importFrom dplyr inner_join
#'
#' @return A list of expression matrices (of different species) with only shared genes
#'
#' @keywords shared gene expression
#'
#' @examples
#'
#' data(NeuroGenesis4)
#' SJDdataIN = sjdWrap(
#' data.list = NeuroGenesis4,
#' species.vector=c("human","human","human","mouse"),
#' geneType.vector=c("symbol","ensembl","symbol","symbol"),
#' geneType.out="symbol",
#' species.out="human")
#'
#' @export

sjdWrap <- function(data.list, species.vector, geneType.vector, geneType.out="symbol", species.out){
    N = length(data.list) #number of datasets

    if(N != length(species.vector)){
        print("length of data.list != length(species.vector)")
        return(NULL)
    }

    if(N != length(geneType.vector)){
        print("length of data.list != length(geneType.vector)")
        return(NULL)
    }


    genes.tbl.list = list()
    new.data.list = list()

    ## Special case that all species and genes are the same
    if(length(unique(species.vector)) == 1 && length(unique(geneType.vector)) == 1){
        gene_names = rownames(data.list[[1]])
        for(i in 1 : N){
            gene_names = intersect(gene_names, rownames(data.list[[i]]))
        }
        for(i in 1 : N){
            new.data.list[[i]] = as.data.frame(data.list[[i]][gene_names, ])
            rownames(new.data.list[[i]]) = gene_names
        }
        names(new.data.list) = names(data.list)
        return(new.data.list)
    }


    cat('Using biomaRt to connect gene IDs across', N, 'datasets:\n')

    for (i in 1 : N) {

        cat('Getting biomaRt IDs for dataset',i,'\n')

        species.1 = species.vector[i]

        genes = rownames(data.list[[i]])

        genes.tbl = getMatch(genes, inSpecies=species.1, inType=geneType.vector[i], newSpecies = species.out)
        genes.tbl = na.omit(genes.tbl)

        genes.tbl$ORIGid = genes.tbl[, "genes"]

        if(geneType.out == "symbol") {
            genes.tbl[, "genes"] = genes.tbl[, 5]
            }#geneType.out

        if(geneType.out == "ensembl"){
            genes.tbl[, "genes"] = genes.tbl[, 6]
            }#geneType.out

        indxNoDup = !duplicated(genes.tbl[, "genes"]) & !is.na(genes.tbl[, "genes"]) & genes.tbl[, "genes"] != ""
        genes.tbl.list[[i]] = genes.tbl[indxNoDup, c("genes", "ORIGid")]
        rm(genes.tbl)
    }

    cat('constructed', length(genes.tbl.list), 'tables of cross-species matching genes\n')

    genes.tbl = inner_join(genes.tbl.list[[1]], genes.tbl.list[[2]], by = 'genes', na_matches = "never") #geneType.out

    if(N >= 3){
        for(i in 3 : N){
            genes.tbl = inner_join(genes.tbl, genes.tbl.list[[i]], by = 'genes', na_matches = "never") #geneType.out
        }
    }

    cat('we found', dim(genes.tbl)[1], 'shared genes in', length(genes.tbl.list), 'datasets\n')

    for (i in 1 : N) {
        indxR = match(genes.tbl[, i + 1], rownames(data.list[[i]])) #ORIGid
        data = data.list[[i]][indxR, ]
        rownames(data) = genes.tbl[, "genes"]
        new.data.list[[i]] = as.data.frame(data)
        rm(data)
    }

    cat('new data list of', length(new.data.list), 'datasets constructed\n')

    names(new.data.list) = names(data.list)

    return(new.data.list)
}
