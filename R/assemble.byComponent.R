#' Assemble Files Based on Component
#'
#' Assemble individual same-component figures from multiple datasets analyzed and plotted by SJD for cross comparison
#'
#' @param SJDScorePlotter.obj A list outputted by the SJDScorePlotter function
#' @param component numer/order of component of interest to print out, i.e. 1 or c(1,2)
#' @param SJD_algorithm SJD_algorithm name of SJD algorithm, i.e. concatICA
#' @param group group name of the weights group from the SJD_algorithm output, i.e 'Shared.All.13'
#'
#' @importFrom stringr str_detect
#'
#' @return A list of images filtered by component
#'
#' @keywords component, images
#'
#' @examples
#'
#' library(ggplot2)
#'
#' data(NeuroGenesis4.afterWrap)
#' data(NeuroGenesis4.info)
#'
#' SampleMetaNamesTable = data.frame(
#'    row.names = names(NeuroGenesis4),
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
#' assemble.byComponent.obj = assemble.byComponent(
#' SJDScorePlotter.obj = SJDScorePlotter.obj,
#' component = c(1, 2),
#' SJD_algorithm = "twoStageLCA",
#' group = 'Shared.All.4')
#'
#' @export


assemble.byComponent <- function(SJDScorePlotter.obj, component, SJD_algorithm, group = NA){
    out_obj_list = list()

    for(comp in component){
        ## Extract names of SJDScorePlotter.obj
        names_plotter_obj = names(SJDScorePlotter.obj)

        names_plotter_obj = names_plotter_obj[which(str_detect(names_plotter_obj, paste0("comp", comp, "of")))]

        ## Filter name based on SJD_algorithm
        names_plotter_obj = names_plotter_obj[which(str_detect(names_plotter_obj, SJD_algorithm))]

        ## Filter name based on group
        if(!is.na(group)){
            names_plotter_obj = names_plotter_obj[which(str_detect(names_plotter_obj, group))]
        }

        for(name in names_plotter_obj){
            out_obj_list[[name]] = SJDScorePlotter.obj[[name]]
        }

    }

    return(out_obj_list)
}

