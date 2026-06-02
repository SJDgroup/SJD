#' Assemble files based on Dataset
#'
#' Assemble individual figures for a single dataset analyzed and plotted by SJD
#'
#' @param SJDScorePlotter.obj A list outputted by the SJDScorePlotter function
#' @param dataset_name dataset/study analyzed by SJD
#' @param SJD_algorithm SJD_algorithm name of SJD algorithm, i.e. concatICA
#' @param group group name of the weights group from the SJD_algorithm output, i.e 'Shared.All.13'
#'
#' @importFrom stringr str_detect
#'
#' @return a list of images filtered by dataset
#'
#' @keywords dataset, images
#'
#' @examples
#'
#' library(ggplot2)
#'
#' data(NeuroGenesis4.afterWrap)
#' data(NeuroGenesis4.info)
#'
#' SampleMetaNamesTable = data.frame(
#'    row.names = names(NeuroGenesis4.afterWrap),
#'    Type = c('Yaxis','Yaxis','2Dscatter','2Dscatter'),
#'    XaxisColumn = c("X","DAYx","tSNE_1","tsne1:ch1"),
#'    YaxisColumn = c("PJDscores","PJDscores","tSNE_2","tsne2:ch1"),
#'    COLaxisColumn = c("color","colorBYlabelsX","PJDscores","PJDscores"),
#'    PCHColumn = c("","","","")
#' )
#'
#' grp = list(
#' Shared.All.4 = c(1 : 4),
#' Shared.bulk.2 = c(1, 2),
#' Shared.sc.2 = c(3, 4),
#' Hs.Meisnr.1 = c(1),
#' Hs.AZ.1 = c(2),
#' Gesch.1 = c(3),
#' Telley.1 = c(4)
#' )
#' dims = c(2, 2, 2, 2, 2, 2, 2)
#'
#' lbb = "NeuroGenesis4.p2"
#'
#' twoStageLCA.out = twoStageLCA(dataset = NeuroGenesis4.afterWrap, group = grp, comp_num = dims)
#'
#' SJDScorePlotter.obj = SJDScorePlotter(
#'     SJDalg = "twoStageLCA",
#'     scores = twoStageLCA.out$score_list,
#'     lbb = lbb,
#'     info = NeuroGenesis4.info,
#'     SampleMetaNamesTable = SampleMetaNamesTable
#' )
#'
#' assemble.byDataset.obj = assemble.byDataset(
#' SJDScorePlotter.obj = SJDScorePlotter.obj,
#' dataset_name = "Meissner.inVitro.bulk.Hs",
#' SJD_algorithm = "twoStageLCA",
#' group = NA)
#'
#' @export

assemble.byDataset <- function(SJDScorePlotter.obj, dataset_name, SJD_algorithm, group = NA){
    ## Extract names of SJDScorePlotter.obj
    names_plotter_obj = names(SJDScorePlotter.obj)

    ## Filter name based on dataset_name
    names_plotter_obj = names_plotter_obj[which(str_detect(names_plotter_obj, paste0("\\.data_", dataset_name, "\\.")))]

    ## Filter name based on SJD_algorithm
    names_plotter_obj = names_plotter_obj[which(str_detect(names_plotter_obj, SJD_algorithm))]

    ## Filter name based on group
    if(!is.na(group)){
        names_plotter_obj = names_plotter_obj[which(str_detect(names_plotter_obj, group))]
    }

    ## Filtered list of obj
    out_obj_list = list()

    for(name in names_plotter_obj){
        out_obj_list[[name]] = SJDScorePlotter.obj[[name]]
    }

    return(out_obj_list)
}
