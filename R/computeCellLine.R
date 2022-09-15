#' @export
#' @importFrom foreach foreach %do%
#' @importFrom rhdf5 h5read
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom octad.db get_ExperimentHub_data
#' @importFrom utils write.csv data txtProgressBar read.csv2 head read.csv
#' @importFrom grDevices pdf
#' @importFrom utils write.csv data txtProgressBar read.csv2 head read.csv
#' @importFrom grDevices pdf
####### computeCellLine #######
computeCellLine <- function(case_id = case_id, expSet = NULL, LINCS_overlaps = TRUE,
                            source = c("octad.small", "octad.whole", "expSet"), file = NULL, output = TRUE,
                            outputFolder = NULL) {
  # STOPS
  if (missing(case_id)) {
    stop("Case ids vector input not found")
  }

  #output check
  if (output==TRUE&is.null(outputFolder)) {
    outputFolder <- tempdir()
    message('outputFolder is NULL, writing output to tempdir()')
  }else if (output==TRUE&!is.null(outputFolder)){
    if (output==TRUE&!dir.exists(outputFolder)) {
      dir.create(outputFolder)
    } else if (output==TRUE&dir.exists(outputFolder)){
      warning('Existing directory ', outputFolder, ' found, containtment might be overwritten')
    }
  }
  
  


  # helper function by Ke
  pick.out.cell.line <- function(expr.of.samples, expr.of.cell.lines, marker.gene) {
    marker.gene <- intersect(rownames(expr.of.samples), (marker.gene))
    marker.gene <- intersect(rownames(expr.of.cell.lines), (marker.gene))
    correlation.matrix <- cor(expr.of.samples[marker.gene, ], expr.of.cell.lines[marker.gene, ], method = "spearman")
    correlation.matrix[is.na(correlation.matrix)] <- 0
    cell.line.median.cor <- apply(correlation.matrix, 2, median) %>%
      sort(decreasing = TRUE)
    best.cell.line <- names(cell.line.median.cor)[1]
    cell.line <- setdiff(names(cell.line.median.cor), best.cell.line)
    p.value.vec <- NULL
    for (i in cell.line) {
      v <- correlation.matrix[, i]
      p.value.vec <- c(p.value.vec, wilcox.test(correlation.matrix[, best.cell.line], v, alternative = "greater", paired = TRUE)$p.value)
    }
    names(p.value.vec) <- setdiff(names(cell.line.median.cor), best.cell.line)
    fdr.vec <- p.adjust(p.value.vec, method = "fdr")
    list(cell.line.median.cor = cell.line.median.cor, best.cell.line = best.cell.line, compare.fdr.vec = fdr.vec, correlation.matrix = correlation.matrix)
  }
  
  
  source=match.arg(source)
  if (source == "octad.whole") {
    message("loading whole octad expression data for ", length(c(case_id)), " samples")
    transcripts <- as.character(rhdf5::h5read(file, "meta/transcripts"))
    samples <- as.character(rhdf5::h5read(file, "meta/samples"))
    case_counts <- rhdf5::h5read(file, "data/count", index = list(seq_len(length(transcripts)), which(samples %in% case_id)))
    colnames(case_counts) <- samples[samples %in% case_id]
    rownames(case_counts) <- transcripts
    case_id <- samples[samples %in% case_id]
    H5close()
    case_counts <- cbind(case_counts)
    message("computing correlation between cell lines and selected samples", sep = " ")
  } else if (source == "octad.small") {
    message("loading whole octad expression data for ", length(c(case_id)), " samples", sep = " ")
    octad.LINCS.counts <- suppressMessages(get_ExperimentHub_data("EH7273"))
    case_counts <- octad.LINCS.counts[, c(case_id)]
    message("computing correlation between cell lines and selected samples", sep = " ")
  } else if (source != "octad"|source != "octad.small" & !missing(expSet)) {
    case_counts <- expSet[, case_id]
  } else {
    stop("Expression data not sourced, please, modify expSet option")
  }

  if (LINCS_overlaps == TRUE) {
    CCLE.overlaps <- get_ExperimentHub_data("EH7262")
    CCLE.median <- apply(CCLE.overlaps, 1, median)
    CCLE.log2.read.count.matrix <- get_ExperimentHub_data("EH7261")
    CCLE.expressed.gene <- names(CCLE.median)[CCLE.median > 1]
  } else {
    CCLE.log2.read.count.matrix <- get_ExperimentHub_data("EH7261") # bioconductor addon
    CCLE.median <- apply(CCLE.log2.read.count.matrix, 1, median)
    CCLE.expressed.gene <- names(CCLE.median)[CCLE.median > 1]
  }
   # bioconductor addon
  tmp <- CCLE.log2.read.count.matrix[CCLE.expressed.gene, ]
  tmp.rank <- apply(tmp, 2, rank)
  rank.mean <- apply(tmp.rank, 1, mean)
  rank.sd <- apply(tmp.rank, 1, sd)
  CCLE.rna.seq.marker.gene.1000 <- names(sort(rank.sd, decreasing = TRUE))[seq_len(1000)]
  TCGA.vs.CCLE.polyA.expression.correlation.result <- pick.out.cell.line(case_counts, CCLE.overlaps, CCLE.rna.seq.marker.gene.1000)
  correlation.dataframe <- TCGA.vs.CCLE.polyA.expression.correlation.result$cell.line.median.cor %>%
    as.data.frame()
  colnames(correlation.dataframe) <- "cor"


  topline <- data.frame(medcor = TCGA.vs.CCLE.polyA.expression.correlation.result$cell.line.median.cor) # could also do first

  if (output == TRUE) {
    write.csv(correlation.dataframe, file = file.path(outputFolder, "CellLineCorrelations.csv"))
    return(topline)
  } else {
    return(topline)
  }
}
