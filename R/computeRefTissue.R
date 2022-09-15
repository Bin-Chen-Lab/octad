#' @export
#' @importFrom magrittr %>%
#' @import octad.db stats
#' @importFrom dplyr select arrange mutate desc
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom octad.db get_ExperimentHub_data
#' @importFrom utils write.csv data txtProgressBar read.csv2 head read.csv
#' @importFrom grDevices pdf
computeRefTissue <- function(case_id = NULL, adjacent = FALSE, source = "octad", n_varGenes = 500, method = c("varGenes",'random'), expSet = NULL,
                             control_size = length(case_id), outputFolder = NULL, cor_cutoff = "0", output = TRUE) {
  if (missing(case_id)) {
    stop("Case ids vector input not found")
  }
  if (source == "octad") {
    expSet <- get_ExperimentHub_data("EH7265")
    case_id <- case_id[case_id %in% colnames(expSet)]
    # if we pick adjacent, filter them out
    if (adjacent == TRUE) {
      phenoDF <- get_ExperimentHub_data("EH7274")
      adjacent_ids <- as.vector(subset(phenoDF, phenoDF$sample.type == "adjacent")$sample.id) # bioconductor replace
      normal_id <- as.vector(subset(phenoDF, phenoDF$sample.type == "normal")$sample.id) # bioconductor replace
      normal_id <- c(adjacent_ids, normal_id)
    } else {
      normal_id <- as.vector(subset(phenoDF, phenoDF$sample.type == "normal")$sample.id) # bioconductor replace
    }
  } else if (source != "octad" & missing(expSet)) {
    stop("expSet is not supported")
  } else if (source != "octad") {
    normal_id <- colnames(expSet)[!colnames(expSet) %in% case_id]
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

  normal_id <- normal_id[normal_id %in% colnames(expSet)]
  
  
  method=match.arg(method)
  if (method == "random") {
    GTEXid <- sample(normal_id, size = control_size)
    return(GTEXid)
  } else if (method == "varGenes") {
    expSet_normal <- expSet[, as.vector(normal_id)]
    expSet_case <- expSet[, as.vector(case_id)]
    # varGenes look at the top varying genes (IQR) within normal tissue expression and varies them to the case tissues
    iqr_gene <- apply(expSet_normal, 1, stats::IQR) # get the IQR per gene
    varying_genes <- order(iqr_gene, decreasing = TRUE)[seq_len(min(n_varGenes, length(iqr_gene)))]

    # get the correlation matrix for each normal id and each case id
    normal_dz_cor <- cor(expSet_normal[varying_genes, ], expSet_case[varying_genes, ], method = "spearman")
    normal_dz_cor_each <- apply(normal_dz_cor, 1, median) # getting the median correlation btw each normal tissue to the case overall
    cor <- data.frame(cor = sort(normal_dz_cor_each, decreasing = TRUE))
    sample.id <- row.names(cor)
    normal_dz_cor_eachDF <- cor %>%
      dplyr::mutate(sample.id) %>%
      dplyr::select(sample.id, cor)
    cutoff <- stats::quantile(normal_dz_cor_eachDF$cor, probs = seq(0, 1, 0.05), na.rm = TRUE)[paste0(cor_cutoff, "%")]

    GTEXid_temp <- subset(normal_dz_cor_eachDF, cor >= cutoff)
    GTEXid_temp <- GTEXid_temp[order(GTEXid_temp$cor, decreasing = TRUE), ]
    GTEXid <- GTEXid_temp$sample.id

    GTEXid <- GTEXid[seq_len(min(control_size, length(GTEXid)))]

    if (output == TRUE) {
      if (nchar(outputFolder) > 0) {
        tryCatch(write.csv(normal_dz_cor, file = file.path(outputFolder, "case_normal_corMatrix.csv")), error = function(c) "failed to write case normal cor matrix csv. Try checking if your outputFolder string is correct or exists")
        tryCatch(write.csv(normal_dz_cor_eachDF, row.names = FALSE, paste0(outputFolder, "/case_normal_median_cor.csv")),
          error = function(c) "failed to write case normal median correlation csv. Try checking if your outputFolder string is correct or exists"
        )
      } else {
        tryCatch(write.csv(normal_dz_cor, file = file.path(outputFolder, "case_normal_corMatrix.csv")), error = function(c) "failed to write case normal cor matrix csv. Try checking if your outputFolder string is correct or exists")
        tryCatch(write.csv(normal_dz_cor_eachDF, row.names = FALSE, paste0(outputFolder, "case_normal_median_cor.csv")),
          error = function(c) "failed to write case normal median correlation csv. Try checking if your outputFolder string is correct or exists"
        )
      }
    }
    return(GTEXid)
  }
}
