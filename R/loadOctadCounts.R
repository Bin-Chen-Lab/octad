#' @export
#' @import rhdf5
#### runsRGES #######
loadOctadCounts <- function(sample_vector = NULL, type = "tpm", file = NULL) {
  if (missing(sample_vector)) {
    stop("sample_vector input not found")
  }
  
  if (is.null(file)) {
    stop("h5 file with counts not provided")
  }
  
  transcripts = as.character(rhdf5::h5read(file, "meta/transcripts"))
  samples = as.character(rhdf5::h5read(file, "meta/samples"))
  if (type == "counts") {
    message("loading", length(transcripts), "log2 expression values for", length(sample_vector), "samples", sep = " ")
    exprData = rhdf5::h5read(file, "data/count", index = list(seq_len(length(transcripts)), which(samples %in% sample_vector)))
    rownames(exprData) = transcripts
    colnames(exprData) = samples[samples %in% sample_vector]
    rhdf5::H5close()
    return(exprData)
  } else if (type == "tpm") {
    message("loading", length(transcripts), "TPM expression values for", length(sample_vector), "samples", sep = " ")
    exprData = rhdf5::h5read(file, "data/tpm", index = list(seq_len(length(transcripts)), which(samples %in% sample_vector)))
    rownames(exprData) = transcripts
    colnames(exprData) = samples[samples %in% sample_vector]
    rhdf5::H5close()
    return(exprData)
  } else {
    message("Some parameters are missing, check everything is ok and restart function")
  }
  rhdf5::H5close()
}