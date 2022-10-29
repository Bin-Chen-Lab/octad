#' @export
#' @importFrom data.table dcast.data.table
#' @importFrom htmlwidgets saveWidget
#' @importFrom reshape2 melt
#' @importFrom plotly ggplotly add_annotations layout as_widget
#' @import octad.db
#' @import ggplot2
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom octad.db get_ExperimentHub_data
#' @importFrom utils write.csv data txtProgressBar read.csv2 head read.csv
#' @importFrom grDevices pdf
####### topLineEval #######
topLineEval <- function(topline = NULL, mysRGES = NULL, outputFolder = NULL) {
    if (missing(mysRGES)) {
        stop("sRGES signature input not found")
    }
    if (length(topline) == 0) {
        stop("No cell lines found in the input")
    }

    if (is.null(outputFolder)) {
        outputFolder <- tempdir()
        message("outputFolder is NULL, writing output to tempdir()")
    }


    toplineName <- paste(topline, collapse = "_")

    cell.line.folder <- file.path(outputFolder, "CellLineEval")
    if (!dir.exists(cell.line.folder)) {
        # check for overwrite
        dir.create(cell.line.folder)
    } else if (dir.exists(cell.line.folder)) {
        # check for overwrite
        warning("Existing directory", cell.line.folder, "found, containtment might be overwritten")
    }




    mysRGES$pert_iname <- toupper(mysRGES$pert_iname)

    CTRPv2.ic50 <- data.table::dcast.data.table(octad.db::get_ExperimentHub_data("EH7264"), drugid ~ cellid, value.var = "ic50_recomputed",
        fun.aggregate = median)  # bioconductor replace

    colnames(CTRPv2.ic50) <- gsub("[^0-9A-Za-z///' ]", "", colnames(CTRPv2.ic50))
    colnames(CTRPv2.ic50) <- toupper(colnames(CTRPv2.ic50))
    colnames(CTRPv2.ic50) <- gsub(" ", "", colnames(CTRPv2.ic50))
    ###
    CTRPv2.ic50 <- as.data.frame(CTRPv2.ic50)
    CTRP.IC50 <- subset(CTRPv2.ic50, select = c("DRUGID", topline))


    CTRP.IC50.m <- as.data.frame(reshape2::melt(CTRP.IC50, id.vars = "DRUGID"))
    CTRP.IC50.medianIC50 <- aggregate(CTRP.IC50.m[3], by = list(CTRP.IC50.m$DRUGID), FUN = median)
    CTRP.IC50.medianIC50$medIC50 <- CTRP.IC50.medianIC50$value
    CTRP.IC50.medianIC50$value <- NULL
    CTRP.IC50.medianIC50$DRUGID <- CTRP.IC50.medianIC50$Group.1
    CTRP.IC50.medianIC50$Group.1 <- NULL
    ###
    CTRP.IC50.medianIC50 <- CTRP.IC50.medianIC50[is.finite(CTRP.IC50.medianIC50$medIC50), ]


    CTRPv2.auc <- data.table::dcast.data.table(get_ExperimentHub_data("EH7264"), drugid ~ cellid, value.var = "auc_recomputed",
        fun.aggregate = median)

    colnames(CTRPv2.auc) <- gsub("[^0-9A-Za-z///' ]", "", colnames(CTRPv2.auc))
    colnames(CTRPv2.auc) <- toupper(colnames(CTRPv2.auc))
    colnames(CTRPv2.auc) <- gsub(" ", "", colnames(CTRPv2.auc))
    ###
    CTRPv2.auc <- as.data.frame(CTRPv2.auc)
    CTRP.auc <- subset(CTRPv2.auc, select = c("DRUGID", topline))
    ###

    ###
    CTRP.auc.m <- as.data.frame(reshape2::melt(CTRP.auc, id.vars = "DRUGID"))
    CTRP.auc.medianauc <- aggregate(CTRP.auc.m[3], by = list(CTRP.auc.m$DRUGID), FUN = median)
    CTRP.auc.medianauc$medauc <- CTRP.auc.medianauc$value
    CTRP.auc.medianauc$value <- NULL
    CTRP.auc.medianauc$DRUGID <- CTRP.auc.medianauc$Group.1
    CTRP.auc.medianauc$Group.1 <- NULL
    ###
    CTRP.auc.medianauc <- CTRP.auc.medianauc[is.finite(CTRP.auc.medianauc$medauc), ]
    mysRGES$pert_iname <- toupper(mysRGES$pert_iname)
    CTRP.IC50.medianIC50$DRUGID <- toupper(CTRP.IC50.medianIC50$DRUGID)
    testdf <- merge(mysRGES, CTRP.IC50.medianIC50, by.x = "pert_iname", by.y = "DRUGID")
    IC50.cortest <- cor.test(testdf$sRGES, log10(testdf$medIC50))
    ic50pval <- IC50.cortest$p.value
    ic50rho <- IC50.cortest$estimate
    mylabel <- c(`p-value` = ic50pval, Rho = ic50rho)

    testdf$StronglyPredicted <- NA
    testdf$StronglyPredicted <- ifelse(testdf$sRGES < -0.2, "Yes", "No")

    StronglyPredicted <- testdf$StronglyPredicted

    Legend.title <- "Strongly <br>Predicted"
    Legend.label1 <- "No"
    Legend.label2 <- "Yes"
    Title <- "Top Line recomputed log(ic50) vs sRGES"
    xaxis <- "sRGES"
    yaxis <- "log(ic50)"
    options(warn = -1)
    p <- ggplot2::ggplot(data = testdf, aes(x = testdf$sRGES, y = log10(testdf$medIC50))) + ggplot2::geom_point(aes(color = StronglyPredicted,
        text = paste("Drug: ", testdf$pert_iname, "<br>sRGES: ", testdf$sRGES))) + ggplot2::geom_smooth(aes(label = ic50rho),
        method = "lm", se = FALSE, color = "black", size = 0.5) + ggplot2::scale_color_discrete(name = Legend.title) +
        ggplot2::labs(x = xaxis, y = yaxis, title = Title) + ggplot2::theme(legend.position = "right", legend.background = element_rect(fill = "#F5F5F5"),
        legend.title = element_blank())

    p1 <- plotly::ggplotly(p, tooltip = c("text", "label")) %>%
        plotly::layout(margin = list(l = 15))  # change to white background

    ic50graph <- p1 %>%
        plotly::add_annotations(text = Legend.title, xref = "paper", yref = "paper", x = 1.02, xanchor = "left",
            y = 0.8, yanchor = "bottom", legendtitle = TRUE, showarrow = FALSE) %>%
        plotly::layout(legend = list(y = 0.8, yanchor = "top"))



    file_name <- file.path(cell.line.folder, paste0(topline, "_ic50_insilico_validation.html"))



    htmlwidgets::saveWidget(plotly::as_widget(ic50graph), file.path(normalizePath(dirname(file_name)), basename(file_name)),
        selfcontained = FALSE)
    # note: files are simply too large to set selfcontained = T. This just causes issues on linux machines.
    # AUC
    testdf2 <- merge(mysRGES, CTRP.auc.medianauc, by.x = "pert_iname", by.y = "DRUGID")
    AUC.cortest <- cor.test(testdf2$sRGES, testdf2$medauc)

    testdf2$StronglyPredicted <- NA
    testdf2$StronglyPredicted <- ifelse(testdf2$sRGES < -0.2, "Yes", "No")
    StronglyPredicted <- testdf2$StronglyPredicted
    aucpval <- AUC.cortest$p.value
    aucrho <- AUC.cortest$estimate

    Legend.title <- "Strongly <br>Predicted"
    Legend.label1 <- "No"
    Legend.label2 <- "Yes"
    Title <- "Top Line recomputed AUC vs sRGES"
    xaxis <- "sRGES"
    yaxis <- "AUC"


    p <- ggplot(testdf2, aes(x = sRGES, y = medauc)) + ggplot2::geom_point(aes(color = StronglyPredicted, text = paste("Drug: ",
        testdf2$pert_iname, "<br>sRGES: ", testdf2$sRGES))) + ggplot2::geom_smooth(aes(label = aucrho), method = "lm",
        se = FALSE, color = "black", size = 0.5) + ggplot2::scale_color_discrete(name = Legend.title) + ggplot2::labs(x = xaxis,
        y = yaxis, title = Title) + ggplot2::theme(legend.position = "right", legend.background = element_rect(fill = "#F5F5F5"),
        legend.title = element_blank())

    p1 <- plotly::ggplotly(p, tooltip = c("text", "label")) %>%
        plotly::layout(margin = list(l = 15))

    aucgraph <- p1 %>%
        plotly::add_annotations(text = Legend.title, xref = "paper", yref = "paper", x = 1.02, xanchor = "left",
            y = 0.8, yanchor = "bottom", legendtitle = TRUE, showarrow = FALSE) %>%
        plotly::layout(legend = list(y = 0.8, yanchor = "top"))

    file_name <- file.path(cell.line.folder, paste0(topline, "_auc_insilico_validation.html"))


    htmlwidgets::saveWidget(plotly::as_widget(aucgraph), file.path(normalizePath(dirname(file_name)), basename(file_name)),
        selfcontained = FALSE)
    # note: files are simply too large to set selfcontained = T. This just causes issues on linux machines.


    # logging cortests
    con <- file(paste0(cell.line.folder, toplineName, "_drug_sensitivity_insilico_results.txt"))
    sink(con, append = TRUE)
    sink(con, append = TRUE, type = "message")

    # print('AUC cortest') print(AUC.cortest) print('IC50 cortest') print(IC50.cortest)

    sink()
    sink(type = "message")
    options(warn = 0)
}
