#' Plot SJD score
#'
#' plot dimensionality reduction scores for each SJD algorithm for dataset analyzed by SJD
#'
#' @param SJDalg SJD algorithm to plot i.e 'twoStageLCA'
#' @param scores score list of the SJD algorithm i.e twoStageLCA$score_list
#' @param lbb dataset label i.e 'NeuroGenesis4'
#' @param info list of sample meta data matrices
#' @param SampleMetaNamesTable dataframe containing column information of each sample meta data matrices
#' @param clrs2end color scale for result scores from other algorithms. Default: c("plum","purple","blue","blue4","black","darkred","red","orange","yellow")
#' @param clrs1end color scale for result scores from sepNMF, concatNMF and jointNMF algorithms. Default: c("black","black","black","darkred","red","orange","yellow")
#'
#' @import plotrix ggplot2 dplyr
#'
#' @return A list containing ggplot object
#'
#' @keywords images
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
#'    PCHColumn = c("","","",""),
#'    inset = c(TRUE, TRUE, TRUE, TRUE),
#'    insetLOC = c("topright", "topright", "topright", "topright"),
#'    insetZoom = c(0.3, 0.3, 0.3, 0.3),
#'    ordDECREASE=c(FALSE, FALSE, FALSE, FALSE),
#'    CLRfoldPRB=c(0.5, 0.5, 0.5, 0.5)
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
#'
#' dims = c(2, 2, 2, 2, 2, 2, 2)
#'
#' twoStageLCA.out = twoStageLCA(dataset = NeuroGenesis4.afterWrap, group = grp, comp_num = dims)
#'
#' SJDScorePlotter.obj = SJDScorePlotter(
#'     SJDalg = "twoStageLCA",
#'     scores = twoStageLCA.out$score_list,
#'     lbb = "NeuroGenesis4.p2",
#'     info = NeuroGenesis4.info,
#'     SampleMetaNamesTable = SampleMetaNamesTable
#' )
#'
#' @export

SJDScorePlotter <- function(
    SJDalg,
    scores,
    lbb,
    info,
    SampleMetaNamesTable,
    clrs2end = c("plum","purple","blue","blue4","black","darkred","red","orange","yellow"),
    clrs1end = c("black","black","black","darkred","red","orange","yellow")
){

    ######################################
    ##
    ## Precheck the input parameters
    ##
    #######################################


    ## Start function

    print("**********************************************************")
    print(paste("plotting scores for ",SJDalg," anaysis.",sep=""))

    if(!SJDalg%in%c("sepNMF","concatNMF","jointNMF","sepPCA","concatPCA","jointPCA","sepICA","concatICA","jointICA","twoStageLCA","twoStageiLCA")){
        print("SJDalg must be one of:")
        print("sepNMF, concatNMF, jointNMF, sepPCA, concatPCA, jointPCA, sepICA, concatICA, jointICA, twoStageLCA, twoStageiLCA")
        return()
    }

    if(SJDalg %in% c("sepNMF","sepPCA","sepICA")){
        SJDsep = TRUE
    }

    if(SJDalg %in% c("concatNMF","jointNMF","concatPCA","jointPCA","concatICA","jointICA","twoStageLCA","twoStageiLCA")){
        SJDsep = FALSE
    }

    if(!identical(names(info), names(scores))){
        print("dataset names in info and scores are not equivalent (or are not in the same order)")
        print("names(info):")
        print(names(info))
        print("names(scores):")
        print(names(scores))
        return()
    }

    ################################
    ##
    ## Compute and plot images
    ##
    ###############################

    dataset_names = names(info)

    ## List to output all images
    image_out_list = list()

    # "separate" SJD analyses have differently structured output from other SJD algorithms:
    if(SJDsep){# BEGIN if(SJDsep)loop

        for(i in 1 : length(dataset_names)){

            dataset_name=dataset_names[i]
            print("*******************************")
            print(paste("dataset : ",dataset_name,sep=""))

            # sepICA, sepPCA，sepNMF different structures happens below here
            score_dimension = dim(scores[[dataset_name]])[1] # dim(sepPCA$score_list[['Hs.AZ']])=2

            if(is.null(SampleMetaNamesTable[dataset_name, 'cexx']) || SampleMetaNamesTable[dataset_name, 'cexx'] == "" || is.na(SampleMetaNamesTable[dataset_name, 'cexx'])){
                cexx = 1
            }else{
                # cexx = strtoi(SampleMetaNamesTable[dataset_name, 'cexx'])
                cexx = as.numeric(SampleMetaNamesTable[dataset_name, 'cexx'])
            }

            ### PREPARE SJD SCORE
            for(k in 1 : score_dimension){#rank loop
                print(k)
                SJDscores = scores[[dataset_name]][k,]
                cpntNM = rownames(scores[[dataset_name]])[k]

                ### PREPARE X Y AXIS
                if(SampleMetaNamesTable[dataset_name,"Type"] == "Yaxis")
                {
                    if(!is.null(SampleMetaNamesTable[dataset_name, "ordDECREASE"])){
                        ord = order(SJDscores,decreasing=SampleMetaNamesTable[dataset_name, "ordDECREASE"])
                    }else{
                        ord = order(SJDscores,decreasing = FALSE)
                    }

                    if(SampleMetaNamesTable[dataset_name,"PCHColumn"] == "" || is.na(SampleMetaNamesTable[dataset_name,"PCHColumn"])){
                        pchh = 19
                    }
                    if(SampleMetaNamesTable[dataset_name,"PCHColumn"] != "" & !is.na(SampleMetaNamesTable[dataset_name, "PCHColumn"])){
                        if(length(unique(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"PCHColumn"])]))<=5 & length(unique(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"PCHColumn"])]))!=1){pchh0=info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"PCHColumn"])];pchh=c(21:25)[as.numeric(as.factor(pchh0))]}
                        if(length(unique(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"PCHColumn"])]))>5 || length(unique(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"PCHColumn"])]))==1){pchh=19}
                    }

                    ## generate image
                    image_out_name = paste("SJDout_", lbb, ".SJDalg_", SJDalg, ".data_", dataset_name, ".comp", k, "of", dim(scores[[dataset_name]])[1], sep="")

                    image_out_list[[image_out_name]] = data.frame(
                        x_axis_value = as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])]),
                        y_axis_value = as.numeric(SJDscores),
                        point_clr = info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"COLaxisColumn"])]
                    ) %>%
                        ggplot(aes(x = x_axis_value, y = y_axis_value)) +
                        geom_point(color = info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"COLaxisColumn"])], cex = cexx, pch = pchh) +
                        xlab(as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])) +
                        ylab(paste0(dataset_name, ",", rownames(scores[[dataset_name]])[k])) +
                        theme_bw() +
                        theme(axis.text=element_text(size = 15),
                              axis.title=element_text(size = 15, face="bold")
                        )

                    if(!is.null(SampleMetaNamesTable[dataset_name, "inset"]) && SampleMetaNamesTable[dataset_name, "inset"] != "" && !is.na(SampleMetaNamesTable[dataset_name, "inset"]) &&  SampleMetaNamesTable[dataset_name, "inset"] == TRUE){
                        insetLOC = SampleMetaNamesTable[dataset_name, "insetLOC"]
                        insetZoom = SampleMetaNamesTable[dataset_name, "insetZoom"]

                        embed_fig = data.frame(SJDscores = SJDscores) %>%
                            ggplot(aes(x = SJDscores)) +
                            geom_density()

                        if(insetLOC == "topleft"){
                            Xmin0 = min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                            Xmax0 = insetZoom * max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord])) + (1 - insetZoom) *min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                            Ymin0 = (1 - insetZoom) * max( as.numeric(SJDscores)) + insetZoom * min(as.numeric(SJDscores));
                            Ymax0 = max(as.numeric(SJDscores));
                        }
                        if(insetLOC == "topright"){
                            Xmin0 = (1 - insetZoom) * max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord])) + insetZoom * min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                            Xmax0 = max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                            Ymin0 = (1 - insetZoom) * max( as.numeric(SJDscores)) + insetZoom * min(as.numeric(SJDscores));
                            Ymax0 = max(as.numeric(SJDscores));
                        }
                        if(insetLOC == "bottomleft"){
                            Xmin0 = min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                            Xmax0 = (1 - insetZoom) * min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord])) + insetZoom * max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                            Ymin0 = min( as.numeric(SJDscores));
                            Ymax0 = (1 - insetZoom) * min( as.numeric(SJDscores)) + insetZoom * max(as.numeric(SJDscores));
                        }
                        if(insetLOC == "bottomright"){
                            Xmin0 = insetZoom * min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord])) + (1 - insetZoom) * max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                            Xmax0 = max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                            Ymin0 = min( as.numeric(SJDscores));
                            Ymax0 = (1 - insetZoom) * min(as.numeric(SJDscores)) + insetZoom * max(as.numeric(SJDscores));
                        }

                        image_out_list[[image_out_name]] = image_out_list[[image_out_name]] +
                            annotation_custom(
                                ggplotGrob(embed_fig),
                                xmin = Xmin0,
                                xmax = Xmax0,
                                ymin = Ymin0,
                                ymax = Ymax0
                            )
                    }
                }

                ### PREPARE X Y AXIS
                if(SampleMetaNamesTable[dataset_name,"Type"] == "2Dscatter")
                {
                    if(!is.null(SampleMetaNamesTable[dataset_name, "ordDECREASE"])){
                        ord = order(SJDscores,decreasing=SampleMetaNamesTable[dataset_name, "ordDECREASE"])
                    }else{
                        ord = order(SJDscores,decreasing = FALSE)
                    }
                    if(SJDalg == "sepNMF" || SJDalg == "concatNMF" || SJDalg == "jointNMF"){
                        clr = color.scale(SJDscores,extremes=clrs1end)
                        CLRfold = FALSE
                    }
                    if(SJDalg!="sepNMF" && SJDalg!="concatNMF" && SJDalg != "jointNMF"){
                        clr = color.scale(SJDscores,extremes=clrs2end)
                        CLRfold = TRUE
                    }
                    if(CLRfold){
                        range(SJDscores)#
                        prb = SampleMetaNamesTable[dataset_name, "CLRfoldPRB"]#.5#.7#.5#.8
                        if(!is.null(SampleMetaNamesTable[dataset_name, "CLRfoldPRB"]) && (SampleMetaNamesTable[dataset_name, "CLRfoldPRB"] != 0 || SampleMetaNamesTable[dataset_name, "CLRfoldPRB"] != 1)){
                            mid = quantile(SJDscores,probs=prb)
                            SJDscores00 = SJDscores-mid
                            if(!is.null(SampleMetaNamesTable[dataset_name, "ordDECREASE"])){
                                ord = order(abs(SJDscores00),decreasing=SampleMetaNamesTable[dataset_name, "ordDECREASE"])
                            }else{
                                ord = order(SJDscores00,decreasing = FALSE)
                            }
                        }
                    }

                    ## Edit size of points
                    if(is.null(SampleMetaNamesTable[dataset_name, 'cexx']) || SampleMetaNamesTable[dataset_name, 'cexx'] == "" || is.na(SampleMetaNamesTable[dataset_name, 'cexx'])){
                        if(length(clr) <= 100){
                            cexx = 2
                        }
                        if(length(clr) > 100 & length(clr) <= 1000){
                            cexx = 1
                        }
                        if(length(clr) > 1000 & length(clr) <= 10000){
                            cexx = .75
                        }
                        if(length(clr) > 10000){
                            cexx = .5
                        }
                    }else{
                        # cexx = strtoi(SampleMetaNamesTable[dataset_name, 'cexx'])
                        cexx = as.numeric(SampleMetaNamesTable[dataset_name, 'cexx'])
                    }
                    if(SampleMetaNamesTable[dataset_name,"PCHColumn"]=="" || is.na(SampleMetaNamesTable[dataset_name,"PCHColumn"])){
                        pchh = 19
                    }
                    if(SampleMetaNamesTable[dataset_name,"PCHColumn"]!="" & !is.na(SampleMetaNamesTable[dataset_name,"PCHColumn"])){
                        if(length(unique(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"PCHColumn"])]))<=5 & length(unique(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"PCHColumn"])])) != 1){
                            pchh0 = info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"PCHColumn"])];
                            pchh = c(21:25)[as.numeric(as.factor(pchh0))][ord]
                        }
                        if(length(unique(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"PCHColumn"])]))>5 || length(unique(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"PCHColumn"])])) == 1){
                            pchh = 19
                        }
                    }

                    ## Save plotted images
                    image_out_name = paste("SJDout_",lbb,".SJDalg_",SJDalg,".data_",dataset_name,".comp",k,"of",dim(scores[[dataset_name]])[1],sep="")

                    image_out_list[[image_out_name]] = data.frame(
                        x_axis_value = as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]),
                        y_axis_value = as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord]),
                        point_clr = clr[ord]
                    ) %>%
                        ggplot(aes(x = x_axis_value, y = y_axis_value)) +
                        geom_point(color = clr[ord], cex = cexx, pch = pchh) +
                        xlab(as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])) +
                        ylab(as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])) +
                        # xlab(rownames(scores[[dataset_name]])[1]) +
                        # ylab(rownames(scores[[dataset_name]])[2]) +
                        ggtitle(paste0(dataset_name, "\n",  rownames(scores[[dataset_name]])[k])) +
                        theme_bw() +
                        theme(axis.text=element_text(size = 15),
                              axis.title=element_text(size = 15, face="bold"))

                    if(!is.null(SampleMetaNamesTable[dataset_name, "inset"]) && SampleMetaNamesTable[dataset_name, "inset"] != "" && !is.na(SampleMetaNamesTable[dataset_name, "inset"]) &&  SampleMetaNamesTable[dataset_name, "inset"] == TRUE){
                        insetLOC = SampleMetaNamesTable[dataset_name, "insetLOC"]
                        insetZoom = SampleMetaNamesTable[dataset_name, "insetZoom"]

                        embed_fig = data.frame(SJDscores = SJDscores) %>%
                            ggplot(aes(x = SJDscores)) +
                            geom_density()

                        if(insetLOC == "topleft"){
                            Xmin0 = min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                            Xmax0 = insetZoom * max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord])) + (1 - insetZoom) *min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                            Ymin0 = (1 - insetZoom) * max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord])) + insetZoom * min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord]));
                            Ymax0 = max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord]));
                        }
                        if(insetLOC == "topright"){
                            Xmin0 = (1 - insetZoom) * max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord])) + insetZoom * min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                            Xmax0 = max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                            Ymin0 = (1 - insetZoom) * max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord])) + insetZoom * min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord]));
                            Ymax0 = max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord]));
                        }
                        if(insetLOC == "bottomleft"){
                            Xmin0 = min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                            Xmax0 = (1 - insetZoom) * min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord])) + insetZoom * max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                            Ymin0 = min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord]));
                            Ymax0 = (1 - insetZoom) * min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord])) + insetZoom * max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord]));
                        }
                        if(insetLOC == "bottomright"){
                            Xmin0 = insetZoom * min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord])) + (1 - insetZoom) * max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                            Xmax0 = max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                            Ymin0 = min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord]));
                            Ymax0 = (1 - insetZoom) * min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord])) + insetZoom * max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord]));
                        }

                        image_out_list[[image_out_name]] = image_out_list[[image_out_name]] +
                            annotation_custom(
                                ggplotGrob(embed_fig),
                                xmin = Xmin0,
                                xmax = Xmax0,
                                ymin = Ymin0,
                                ymax = Ymax0
                            )
                    }

                }
            }
        }
    }# END if(SJDsep)loop

    if(!SJDsep){# BEGIN if(!SJDsep) loop

        for(i in 1 : length(dataset_names)){
            dataset_name = dataset_names[i]
            print("*******************************")
            print(paste("dataset: ",dataset_name, sep=""))

            if(is.null(SampleMetaNamesTable[dataset_name, 'cexx']) || SampleMetaNamesTable[dataset_name, 'cexx'] == "" || is.na(SampleMetaNamesTable[dataset_name, 'cexx'])){
                cexx = 1
            }else{
                # cexx = strtoi(SampleMetaNamesTable[dataset_name, 'cexx'])
                cexx = as.numeric(SampleMetaNamesTable[dataset_name, 'cexx'])
            }

            for (j in 1 : length(scores[[dataset_name]])) {

                print(paste("grouping: ", names(scores[[dataset_name]])[j], sep=""))

                if (is.na(scores[[dataset_name]][[j]][1])) {
                    print(paste("No SJD scores for grouping : ",names(scores[[dataset_name]])[j],", moving to the next SJD output for ",dataset_name,sep=""))
                    next
                }

                score_dimension=dim(scores[[dataset_name]][[j]])[1] # dim(twoStageLCA$score_list[['Hs.AZ']])=2

                ### PREPARE SJD SCORE
                for(k in 1:score_dimension) {#rank loop
                    print(k)
                    SJDscores = scores[[dataset_name]][[j]][k,]
                    cpntNM = rownames(scores[[dataset_name]][[j]])[k]

                    ### PREPARE X Y AXIS
                    if(SampleMetaNamesTable[dataset_name,"Type"] == "Yaxis")
                    {
                        if(!is.null(SampleMetaNamesTable[dataset_name, "ordDECREASE"])){
                            ord = order(SJDscores,decreasing=SampleMetaNamesTable[dataset_name, "ordDECREASE"])
                        }else{
                            ord = order(SJDscores,decreasing = FALSE)
                        }

                        ## change the size of points
                        pchh = 19
                        if(SampleMetaNamesTable[dataset_name,"PCHColumn"] != "" & !is.na(SampleMetaNamesTable[dataset_name,"PCHColumn"])){
                            if(length(unique(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"PCHColumn"])])) <= 5 & length(unique(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"PCHColumn"])])) != 1){
                                pchh0 = info[[dataset_name]][, as.character(SampleMetaNamesTable[dataset_name,"PCHColumn"])]
                                pchh = c(21 : 25)[as.numeric(as.factor(pchh0))]
                            }
                        }

                        ## generate image
                        image_out_name = paste("SJDout_",lbb,".SJDalg_",SJDalg,".grp_",names(scores[[dataset_name]])[j], ".data_", dataset_name, ".comp", k, "of", dim(scores[[dataset_name]][[j]])[1], sep="")

                        image_out_list[[image_out_name]] = data.frame(
                            x_axis_value = as.numeric(info[[dataset_name]][, as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])]),
                            y_axis_value = as.numeric(scores[[dataset_name]][[j]][k, ]),
                            point_clr = info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"COLaxisColumn"])]
                        ) %>%
                            ggplot(aes(x = x_axis_value, y = y_axis_value)) +
                            geom_point(color = info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"COLaxisColumn"])], cex = cexx, pch = pchh) +
                            xlab(as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])) +
                            ylab(paste0(dataset_name, ",", rownames(scores[[dataset_name]][[j]])[k])) +
                            theme_bw() +
                            theme(axis.text=element_text(size = 15),
                                  axis.title=element_text(size = 15, face="bold")
                            )

                        if(!is.null(SampleMetaNamesTable[dataset_name, "inset"]) && SampleMetaNamesTable[dataset_name, "inset"] != "" && !is.na(SampleMetaNamesTable[dataset_name, "inset"]) &&  SampleMetaNamesTable[dataset_name, "inset"] == TRUE){
                            insetLOC = SampleMetaNamesTable[dataset_name, "insetLOC"]
                            insetZoom = SampleMetaNamesTable[dataset_name, "insetZoom"]

                            embed_fig = data.frame(SJDscores = SJDscores) %>%
                                ggplot(aes(x = SJDscores)) +
                                geom_density()

                            if(insetLOC == "topleft"){
                                Xmin0 = min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                                Xmax0 = insetZoom * max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord])) + (1 - insetZoom) *min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                                Ymin0 = (1 - insetZoom) * max( as.numeric(SJDscores)) + insetZoom * min(as.numeric(SJDscores));
                                Ymax0 = max(as.numeric(SJDscores));
                            }
                            if(insetLOC == "topright"){
                                Xmin0 = (1 - insetZoom) * max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord])) + insetZoom * min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                                Xmax0 = max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                                Ymin0 = (1 - insetZoom) * max( as.numeric(SJDscores)) + insetZoom * min(as.numeric(SJDscores));
                                Ymax0 = max(as.numeric(SJDscores));
                            }
                            if(insetLOC == "bottomleft"){
                                Xmin0 = min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                                Xmax0 = (1 - insetZoom) * min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord])) + insetZoom * max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                                Ymin0 = min( as.numeric(SJDscores));
                                Ymax0 = (1 - insetZoom) * min( as.numeric(SJDscores)) + insetZoom * max(as.numeric(SJDscores));
                            }
                            if(insetLOC == "bottomright"){
                                Xmin0 = insetZoom * min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord])) + (1 - insetZoom) * max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                                Xmax0 = max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                                Ymin0 = min( as.numeric(SJDscores));
                                Ymax0 = (1 - insetZoom) * min(as.numeric(SJDscores)) + insetZoom * max(as.numeric(SJDscores));
                            }

                            image_out_list[[image_out_name]] = image_out_list[[image_out_name]] +
                                annotation_custom(
                                    ggplotGrob(embed_fig),
                                    xmin = Xmin0,
                                    xmax = Xmax0,
                                    ymin = Ymin0,
                                    ymax = Ymax0
                                )
                        }
                    }

                    ### PREPARE X Y AXIS
                    if(SampleMetaNamesTable[dataset_name,"Type"] == "2Dscatter")
                    {
                        if(!is.null(SampleMetaNamesTable[dataset_name, "ordDECREASE"])){
                            ord = order(SJDscores,decreasing=SampleMetaNamesTable[dataset_name, "ordDECREASE"])
                        }else{
                            ord = order(SJDscores,decreasing = FALSE)
                        }
                        if(SJDalg == "sepNMF" || SJDalg == "concatNMF" || SJDalg == "jointNMF"){
                            clr = color.scale(SJDscores,extremes=clrs1end)
                            CLRfold = FALSE
                        }

                        if(SJDalg!="sepNMF" & SJDalg != "concatNMF" & SJDalg!="jointNMF")
                        {
                            clr = color.scale(SJDscores,extremes = clrs2end)
                            CLRfold = TRUE
                        }
                        if(CLRfold){
                            range(SJDscores)#
                            prb = SampleMetaNamesTable[dataset_name, "CLRfoldPRB"]#.5#.7#.5#.8
                            if(!is.null(SampleMetaNamesTable[dataset_name, "CLRfoldPRB"]) && (SampleMetaNamesTable[dataset_name, "CLRfoldPRB"]!=0 || SampleMetaNamesTable[dataset_name, "CLRfoldPRB"]!=1)){
                                mid = quantile(SJDscores,probs=prb)
                                SJDscores00 = SJDscores-mid
                                if(!is.null(SampleMetaNamesTable[dataset_name, "ordDECREASE"])){
                                    ord = order(abs(SJDscores00),decreasing=SampleMetaNamesTable[dataset_name, "ordDECREASE"])
                                }else{
                                    ord = order(SJDscores00,decreasing = FALSE)
                                }
                            }
                        }

                        ## Edit size of points
                        if(is.null(SampleMetaNamesTable[dataset_name, 'cexx']) || SampleMetaNamesTable[dataset_name, 'cexx'] == "" || is.na(SampleMetaNamesTable[dataset_name, 'cexx'])){
                            if(length(clr) <= 100){
                                cexx = 2
                            }
                            if(length(clr) > 100 & length(clr) <= 1000){
                                cexx = 1
                            }
                            if(length(clr) > 1000 & length(clr) <= 10000){
                                cexx = .75
                            }
                            if(length(clr) > 10000){
                                cexx = .5
                            }
                        }else{
                            # cexx = strtoi(SampleMetaNamesTable[dataset_name, 'cexx'])
                            cexx = as.numeric(SampleMetaNamesTable[dataset_name, 'cexx'])
                        }

                        ## Edit point type
                        if(SampleMetaNamesTable[dataset_name,"PCHColumn"]=="" || is.na(SampleMetaNamesTable[dataset_name,"PCHColumn"])){
                            pchh = 19
                        }
                        if(SampleMetaNamesTable[dataset_name,"PCHColumn"]!="" & !is.na(SampleMetaNamesTable[dataset_name,"PCHColumn"])){
                            if(length(unique(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"PCHColumn"])]))<=5 & length(unique(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"PCHColumn"])])) != 1){
                                pchh0 = info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"PCHColumn"])]
                                pchh = c(21:25)[as.numeric(as.factor(pchh0))][ord]}

                            if(length(unique(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"PCHColumn"])]))>5 || length(unique(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"PCHColumn"])])) == 1){
                                pchh = 19
                            }
                        }

                        ## Save plotted images

                        image_out_name = paste("SJDout_",lbb,".SJDalg_",SJDalg,".grp_",names(scores[[dataset_name]])[j],".data_",dataset_name,".comp",k,"of",dim(scores[[dataset_name]][[j]])[1],sep="")

                        image_out_list[[image_out_name]] = data.frame(
                            x_axis_value = as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]),
                            y_axis_value = as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord]),
                            point_clr = clr[ord]
                        ) %>%
                            ggplot(aes(x = x_axis_value, y = y_axis_value)) +
                            geom_point(color = clr[ord], cex = cexx, pch = pchh) +
                            xlab(as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])) +
                            ylab(as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])) +
                            # xlab(rownames(scores[[dataset_name]][[j]])[1]) +
                            # ylab(rownames(scores[[dataset_name]][[j]])[2]) +
                            ggtitle(paste0(dataset_name, "\n", rownames(scores[[dataset_name]][[j]])[k])) +
                            theme_bw() +
                            theme(axis.text=element_text(size = 15),
                                  axis.title=element_text(size = 15, face="bold"))#+
                        # inset color legend:
                        # + annotation_custom(ggplotGrob(p2inset), xmin = Xmin0, xmax = Xmax0, ymin = Ymin0, ymax = Ymax0)

                        if(!is.null(SampleMetaNamesTable[dataset_name, "inset"]) && SampleMetaNamesTable[dataset_name, "inset"] != "" && !is.na(SampleMetaNamesTable[dataset_name, "inset"]) &&  SampleMetaNamesTable[dataset_name, "inset"] == TRUE){
                            insetLOC = SampleMetaNamesTable[dataset_name, "insetLOC"]
                            insetZoom = SampleMetaNamesTable[dataset_name, "insetZoom"]

                            embed_fig = data.frame(SJDscores = SJDscores) %>%
                                ggplot(aes(x = SJDscores)) +
                                geom_density()

                            if(insetLOC == "topleft"){
                                Xmin0 = min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                                Xmax0 = insetZoom * max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord])) + (1 - insetZoom) *min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                                Ymin0 = (1 - insetZoom) * max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord])) + insetZoom * min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord]));
                                Ymax0 = max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord]));
                            }
                            if(insetLOC == "topright"){
                                Xmin0 = (1 - insetZoom) * max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord])) + insetZoom * min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                                Xmax0 = max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                                Ymin0 = (1 - insetZoom) * max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord])) + insetZoom * min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord]));
                                Ymax0 = max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord]));
                            }
                            if(insetLOC == "bottomleft"){
                                Xmin0 = min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                                Xmax0 = (1 - insetZoom) * min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord])) + insetZoom * max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                                Ymin0 = min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord]));
                                Ymax0 = (1 - insetZoom) * min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord])) + insetZoom * max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord]));
                            }
                            if(insetLOC == "bottomright"){
                                Xmin0 = insetZoom * min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord])) + (1 - insetZoom) * max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                                Xmax0 = max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"XaxisColumn"])][ord]));
                                Ymin0 = min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord]));
                                Ymax0 = (1 - insetZoom) * min(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord])) + insetZoom * max(as.numeric(info[[dataset_name]][,as.character(SampleMetaNamesTable[dataset_name,"YaxisColumn"])][ord]));
                            }

                            image_out_list[[image_out_name]] = image_out_list[[image_out_name]] +
                                annotation_custom(
                                    ggplotGrob(embed_fig),
                                    xmin = Xmin0,
                                    xmax = Xmax0,
                                    ymin = Ymin0,
                                    ymax = Ymax0
                                )

                            # Xrng=abs(par('usr')[1]-par('usr')[2]);Yrng=abs(par('usr')[3]-par('usr')[4])
                            # p2inset=
                            # plot(density(SJDscores),xlab="",ylab="",main="",xaxt="n",yaxt="n")
                            # plot(x=SJDscores,y=rep(0,length(SJDscores)),col=clr,pch="|",xlab="",ylab="",main="",xaxt="n",yaxt="n",bty="n");par(plt=plt0)
                        }

                    }
                }
            }
        }
        return(image_out_list)
    }# END if(!SJDsep) loop
}

