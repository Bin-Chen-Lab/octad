#' Open Cancer TherApeutic Discovery (OCTAD) database package
#'
#' Open Cancer TherApeutic Discovery (OCTAD) package implies sRGES approach for the drug discovery. 
#' The essential idea is to identify drugs that reverse the gene expression signature of a disease by tamping down over-expressed genes and stimulating weakly expressed ones. 
#' The following package contains all required precomputed data for whole OCTAD pipeline computation. 
#' 
#' @section Details:
#'
#'The main functions are:
#'
#' \itemize{
#' \item \code{\link{computeRefTissue}} - Compute reference control samples from OCTAD database using precomputed \code{EncoderDF} models. 
#' \item \code{\link{diffExp}} -  Compute differential expression for case vs control samples. Will produce the file \code{computedEmpGenes.csv} listing empiricaly differentially expressed genes used for RNA-Seq normalization.
#' \item \code{\link{runsRGES}} - Compute sRGES, a score indicating the reveral potency of each drug. It first computes RGES (Reverse Gene Expression Score) for individual instances and then summarizes RGES of invididual drugs (one drug may have multiple instances under different treatment conditions). 
#' \item \code{\link{computeCellLine}} - Compute Correlation between cell lines and vector of case ids.
#' \item \code{\link{topLineEval}} - Evaluate predictions using pharmacogenomics data. Given a cell line, the function computes the correlation between sRGES and drug sensitivity data taken from CTRP. A higher correlation means a better prediction. The cell line could be computed from computeCellLine.
#' \item \code{\link{octadDrugEnrichment}} - Perform enrichment analysis of drug hits based on chemical structures, drug-targets, and pharmacological classifications. An enrichment score calculated using ssGSEA and a p-value computed through a permutation test are provided. 
#' }
#'
#'For detailed information on usage, see the package vignette, by typing
#' \code{vignette('octad')}, or the workflow linked to on the first page
#' of the vignette.
#' 
#' 
#' The code can be viewed at the GitHub repository,
#' which also lists the contributor code of conduct:
#'
#' \url{https://github.com/Bin-Chen-Lab/OCTAD}
#' 
#' 
#' @references
#' Zeng, B., Glicksberg, B.S., Newbury, P., Chekalin, E., Xing, J., Liu, K., Wen, A., Chow, C. and Chen, B., 2021. OCTAD: an open workspace for virtually screening therapeutics targeting precise cancer patient groups using gene expression features. Nature protocols, 16(2), pp.728-753.
#' [https://www.nature.com/articles/s41596-020-00430-z](https://www.nature.com/articles/s41596-020-00430-z)
#' @docType package
#' @name octad
NULL