#' Species matching function
#'
#' Use biomaRt to retrieve additional gene identifiers for a vector of gene identifiers in one species (e.g. 'human'), or gene identifiers for orthologous genes in a second species (e.g. 'mouse')
#'
#' @param genes character vector of gene identifiers in either ensembl gene identifier or gene symbol format
#' @param inSpecies a character vector indicating the species from which 'genes' come; must be one of: "human", "mouse", "roundworm", "fruitfly", "zebrafish", "chicken", "rat", "guinea pig", "golden hamster", "rabbit","pig", "sheep", "cow", "dog", "cat", "macaque", "bonobo", "chimpanzee"
#' @param inType a single character value indicating the type of gene identifiers being passed in the "genes" argument; must be one of: "symbol", "ensembl"
#' @param newSpecies a character vector indicating the species for which orthologous gene identifiers are desired; must be one of: "human", "mouse", "roundworm", "fruitfly", "zebrafish", "chicken", "rat", "guinea pig", "golden hamster", "rabbit", "pig", "sheep", "cow", "dog", "cat", "macaque", "bonobo", "chimpanzee"
#' @param moreAttrIn character vector of other gene attributes that you want returned for the input species from biomaRt - to add this argument you must know the names of the fields in the species-specific biomaRt that you are requesting. default is NA.
#' @param moreAttrNew character vector of other gene attributes that you want returned for the output species from biomaRt - to add this argument you must know the names of the fields in the species-specific biomaRt that you are requesting. default is NA.
#' @param useNewestVersion logical indicating if the function should attempt to use the latest version of ensembl, or to use the Aug 2020 archive version that is more stable. default is FALSE.
#'
#' @importFrom biomaRt useMart getLDS
#'
#' @return A matrix whose first column contains the exact gene identifiers (and in the same order) that were sent to the function in the "genes" argument, followed by 3 columns of gene identifiers that were retrieved from biomaRt from both the "inSpecies" and the "newSpecies" (for each species, these 3 identifiers are: gene symbol, ensembl gene ID, and text description).
#'
#' @keywords shared genes
#'
#' @examples
#'
#' data(NeuroGenesis4)
#' out = getMatch(
#' rownames(NeuroGenesis4$Meissner.inVitro.bulk.Hs),
#' inSpecies = 'human',
#' inType = 'symbol',
#' newSpecies = 'mouse')
#'
#' @export

getMatch <- function(genes, inSpecies, inType, newSpecies, useNewestVersion = FALSE, moreAttrIn = NA, moreAttrNew = NA){
    species = data.frame(
        species.nm = c("human", "mouse", "roundworm",
                       "fruitfly", "zebrafish", "frog", "chicken", "rat", "guinea pig",
                       "golden hamster", "rabbit", "pig", "sheep", "cow", "dog",
                       "cat", "macaque", "bonobo", "chimpanzee"),
        species.latin = c("homo sapiens",
                          "mus musculus", "caenorhabditis elegans", "drosphila melanagaster",
                          "danio rerio", "xenopus tropicalis", "gallus gallus", "ratus norvegicus", "cavia porcellus",
                          "melanochromis auratus", "oryctolagus cuniculus", "sus scrofa domesticus",
                          "ovis aries", "bos taurus", "canis lupus familiaris",
                          "felis catus", "macaca mulatta", "pan paniscus", "pan troglodytes"),
        ensembl.nms = c("hsapiens_gene_ensembl", "mmusculus_gene_ensembl",
                        "celegans_gene_ensembl", "dmelanogaster_gene_ensembl",
                        "drerio_gene_ensembl", "xtropicalis_gene_ensembl", "ggallus_gene_ensembl", "rnorvegicus_gene_ensembl",
                        "cporcellus_gene_ensembl", "mauratus_gene_ensembl",
                        "ocuniculus_gene_ensembl", "sscrofa_gene_ensembl",
                        "oaries_gene_ensembl", "btaurus_gene_ensembl", "clfamiliaris_gene_ensembl",
                        "fcatus_gene_ensembl", "mmulatta_gene_ensembl", "ppaniscus_gene_ensembl",
                        "ptroglodytes_gene_ensembl"),
        genesymbol.attr = c("hgnc_symbol",
                            "mgi_symbol", "external_gene_name", "external_gene_name",
                            "hgnc_symbol", "hgnc_symbol", "hgnc_symbol", "external_gene_name",
                            "hgnc_symbol", "mgi_symbol", "hgnc_symbol", "hgnc_symbol",
                            "hgnc_symbol", "hgnc_symbol", "hgnc_symbol", "hgnc_symbol",
                            "hgnc_symbol", "hgnc_symbol", "hgnc_symbol"), row.names = 1)
    mart.in = species[inSpecies, ]$ensembl.nms
    mart.new = species[newSpecies, ]$ensembl.nms

    if(useNewestVersion){
        df.in = useMart(biomart = "ensembl", dataset = mart.in)
        df.new = useMart(biomart = "ensembl", dataset = mart.new)
    }else{
        df.in = useMart(biomart = "ensembl", dataset = mart.in, host = "https://aug2020.archive.ensembl.org")
        df.new = useMart(biomart = "ensembl", dataset = mart.new, host = "https://aug2020.archive.ensembl.org")
    }

    symbol.in = as.character(species[inSpecies, "genesymbol.attr"])
    symbol.new = as.character(species[newSpecies, "genesymbol.attr"])
    if (inType == "symbol") {
        filter = as.character(species[inSpecies, "genesymbol.attr"])
    }
    else if (inType == "ensembl") {
        filter = "ensembl_gene_id"
    }
    cat("You have input ", length(genes), " genes\n")
    atb.in = c(symbol.in, "ensembl_gene_id", "description")
    atb.new = c(symbol.new, "ensembl_gene_id", "description")
    if (!is.na(moreAttrIn)) {
        atb.in = c(atb.in, moreAttrIn)
    }
    if (!is.na(moreAttrNew)) {
        atb.new = c(atb.new, moreAttrNew)
    }
    tbl.match = getLDS(attributes = atb.in, mart = df.in, filters = filter,
                       values = genes, martL = df.new, attributesL = atb.new)
    tbl.match[, "Gene.description"] = gsub(" \\[.*", "", tbl.match[,
                                                                   "Gene.description"])
    tbl.match[, "Gene.description.1"] = gsub(" \\[.*", "", tbl.match[,
                                                                     "Gene.description.1"])
    colnames(tbl.match)[1:length(atb.in)] = paste0(colnames(tbl.match)[1:length(atb.in)],
                                                   ".", inSpecies)
    colnames(tbl.match)[(length(atb.in) + 1):(length(atb.in) +
                                                  length(atb.new))] = paste0(colnames(tbl.match)[(length(atb.in) +
                                                                                                      1):(length(atb.in) + length(atb.new))], ".", newSpecies)
    cat("We found ", dim(tbl.match)[1], " matches\n")
    cat(sum(duplicated(tbl.match[, 1])), " of those are duplicates and only keeping the 1st of each\n")
    if (inType == "symbol") {
        tbl.match = tbl.match[match(genes, tbl.match[, 1]), ]
        rownames(tbl.match) = NULL
        tbl.return = cbind(genes, tbl.match, stringsAsFactors = FALSE)
    }
    else if (inType == "ensembl") {
        tbl.match = tbl.match[match(genes, tbl.match[, 2]), ]
        rownames(tbl.match) = NULL
        tbl.return = cbind(genes, tbl.match, stringsAsFactors = FALSE)
    }
    return(tbl.return)
}
