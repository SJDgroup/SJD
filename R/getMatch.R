#' Species matching function
#'
#' Use remart (originally used biomaRt) to retrieve additional gene identifiers for a vector of gene identifiers in one species (e.g. 'human'), or gene identifiers for orthologous genes in a second species (e.g. 'mouse')
#'
#' @param genes character vector of gene identifiers in either ensembl gene identifier or gene symbol format
#' @param inSpecies a character vector indicating the species from which 'genes' come; must be one of: "human", "mouse", "roundworm", "fruitfly", "zebrafish", "chicken", "rat", "guinea pig", "golden hamster", "rabbit","pig", "sheep", "cow", "dog", "cat", "macaque", "bonobo", "chimpanzee"
#' @param inType a single character value indicating the type of gene identifiers being passed in the "genes" argument; must be one of: "symbol", "ensembl"
#' @param newSpecies a character vector indicating the species for which orthologous gene identifiers are desired; must be one of: "human", "mouse", "roundworm", "fruitfly", "zebrafish", "chicken", "rat", "guinea pig", "golden hamster", "rabbit", "pig", "sheep", "cow", "dog", "cat", "macaque", "bonobo", "chimpanzee"
#' @param moreAttrIn character vector of other gene attributes that you want returned for the input species from biomaRt - to add this argument you must know the names of the fields in the species-specific biomaRt that you are requesting. default is NA.
#' @param moreAttrNew character vector of other gene attributes that you want returned for the output species from biomaRt - to add this argument you must know the names of the fields in the species-specific biomaRt that you are requesting. default is NA.
#' @param forceINPUTorder1to1 logical indicating if you want the result to be in the same order + structure as the "genes" input. TRUE by default.
#'
#' @importFrom remart getLDS getBM
#'
#' @return A matrix whose first column contains the exact gene identifiers (and in the same order) that were sent to the function in the "genes" argument, followed by 3 columns of gene identifiers that were retrieved from biomaRt from both the "inSpecies" and the "newSpecies" (for each species, these 3 identifiers are: gene symbol, ensembl gene ID, and text description).
#'
#' @keywords shared genes
#'
#' @examples
#' \dontrun{
#' data(NeuroGenesis4)
#' out = getMatch(
#'   rownames(NeuroGenesis4$Meissner.inVitro.bulk.Hs),
#'   inSpecies = 'human',
#'   inType = 'symbol',
#'   newSpecies = 'mouse'
#' )
#' }
#'
#' @export

getMatch <- function(genes, inSpecies, inType, newSpecies, useNewestVersion = FALSE, moreAttrIn = NA, moreAttrNew = NA, forceINPUTorder1to1=TRUE){

	species = data.frame(
        species.nm = c("human", "mouse", "roundworm",
                       "fruitfly", "zebrafish", "frog", "chicken", "rat", "guinea pig",
                       "golden hamster", "rabbit", "pig", "sheep", "cow", "dog",
                       "cat", "macaque", "bonobo", "chimpanzee","ferret"),
        species.latin = c("homo sapiens",
                          "mus musculus", "caenorhabditis elegans", "drosphila melanagaster",
                          "danio rerio", "xenopus tropicalis", "gallus gallus", "ratus norvegicus", "cavia porcellus",
                          "melanochromis auratus", "oryctolagus cuniculus", "sus scrofa domesticus",
                          "ovis aries", "bos taurus", "canis lupus familiaris",
                          "felis catus", "macaca mulatta", "pan paniscus", "pan troglodytes","Mustela putorius furo")
		, row.names = 1)
		
    if (inType == "symbol") {filter = "external_gene_name"}
    else if (inType == "ensembl") {filter = "ensembl_gene_id"}
	
    cat("You have input ", length(genes), " genes\n")
    atb.in = c("external_gene_name", "ensembl_gene_id", "description")
    atb.new = c("external_gene_name", "ensembl_gene_id", "description")
	
    if (!is.na(moreAttrIn)) {atb.in = c(atb.in, moreAttrIn)}
    if (!is.na(moreAttrNew)) {atb.new = c(atb.new, moreAttrNew)}
    
    Ngs=length(genes)
	Nks=ceiling(Ngs/1000)#number of 1k chunks we need to send becasue Ensembl REST API allows only 1k requests!!!!????
	cat("N genes =", Ngs,", so sending", Nks, "chunks of 1,000 genes to Ensembl REST API via remart (1k is the limit). This may take a while ...\n")
	
	if(inSpecies!=newSpecies){
		if(length(genes)<=1000){tbl.match = getLDS(attributes = atb.in, species = inSpecies, filters = filter, values = genes, speciesL = newSpecies, attributesL = atb.new)}
		if(length(genes)>1000){
			for(ii in 1:Nks){
				if(ii==1){tbl.match = getLDS(attributes = atb.in, species = inSpecies, filters = filter, values = genes[1:1000], speciesL = newSpecies, attributesL = atb.new)}
				if(ii>1 & ii<Nks){tbl.match=rbind(tbl.match,getLDS(attributes = atb.in, species = inSpecies, filters = filter, values = genes[((1000*ii)-999):(1000*ii)], speciesL = newSpecies, attributesL = atb.new))}
				if(ii==Nks){tbl.match=rbind(tbl.match,getLDS(attributes = atb.in, species = inSpecies, filters = filter, values = genes[((1000*ii)-999):Ngs], speciesL = newSpecies, attributesL = atb.new))}
				}
		}
	}
	
	if(inSpecies==newSpecies){
		if(length(genes)<=1000){tbl.match = getBM(attributes = atb.in, species = inSpecies, filters = filter, values = genes)}
		if(length(genes)>1000){
			for(ii in 1:Nks){
				if(ii==1){tbl.match = getBM(attributes = atb.in, species = inSpecies, filters = filter, values = genes[1:1000])}
				if(ii>1 & ii<Nks){tbl.match=rbind(tbl.match,getBM(attributes = atb.in, species = inSpecies, filters = filter, values = genes[((1000*ii)-999):(1000*ii)]))}
				if(ii==Nks){tbl.match=rbind(tbl.match,getBM(attributes = atb.in, species = inSpecies, filters = filter, values = genes[((1000*ii)-999):Ngs]))}
				}
		}
	}

    cat("We found ", dim(tbl.match)[1], " matches\n")
	
    if(forceINPUTorder1to1){
		if(inType=="symbol"){indx=match(genes,tbl.match[,"external_gene_name"])}
		if(inType=="ensembl"){indx=match(genes,tbl.match[,"ensembl_gene_id"])}
		tbl.match=tbl.match[indx,]
	}
	if(!forceINPUTorder1to1){cat("!forceINPUTorder1to1 - not yet implemented - this could be useful for many to one mappings such as fish to mammal.\n")}

	tbl.match[, "description"] = gsub(" \\[.*", "", tbl.match[, "description"])
    if(inSpecies!=newSpecies){
		tbl.match[, "description.1"] = gsub(" \\[.*", "", tbl.match[, "description.1"])
		colnames(tbl.match)[1:length(atb.in)] = paste0(colnames(tbl.match)[1:length(atb.in)], ".", inSpecies)
    	colnames(tbl.match)[(length(atb.in) + 1):(length(atb.in) + length(atb.new))] = paste0(colnames(tbl.match)[(length(atb.in) + 1):(length(atb.in) + length(atb.new))], ".", newSpecies)
	}

	##### have to do some more work and link to "!forceINPUTorder1to1" above to say this:
    #cat(sum(duplicated(tbl.match[, 1])), " of those are duplicates and only keeping the 1st of each\n")
	
    rownames(tbl.match) = NULL
    tbl.match = cbind(genes, tbl.match, stringsAsFactors = FALSE)
	
    return(tbl.match)
}
