\name{runsRGES}
\alias{runsRGES}
\title{Compute sRGES}
\description{
  Compute sRGES, a score indicating the reveral potency of each drug. It first computes RGES (Reverse Gene Expression Score) for individual instances and then summarizes RGES of invididual drugs (one drug may have multiple instances under different treatment conditions). The function returns RGES data.frame containing scores and p.values for every instance. The function also returns \code{data.frame} with drug id in \code{pert_iname} collumn, \code{n}  contains the number of instances for this drug, \code{mean}, \code{median} and \code{sd} of  \code{sRGES} RGES sores. Besides, a number of additional files in the sourced directory: \code{dz_sig_used.csv} contains genes in the disease signature used for computing reverse gene expression scores, \code{sRGES.csv} contains the same data as returned data.frame and \code{all__lincs_score.csv} includes information of RGES.
}
\usage{
 runsRGES(dz_signature,max_gene_size=50, 
          permutations=10000)
}
\arguments{
  \item{dz_signature}{disease signature. Make sure input data frame has a gene \code{Symbol} column, otherwise an error is produced. It must be an UPPERCASE gene symbol.}
\item{choose_fda_drugs}{if \code{TRUE}, only FDA approved drugs are used.}
\item{max_gene_size}{maximum number of disease genes used for drug prediction. By default 50 for each side (up/down).}
\item{cells}{cell ids in \code{lincs_sig_info} file used for prediction. By default, all cell lines are used.}
\item{outputFolder}{folder path to store drug results, by default write results to working directory.}
\item{weight_cell_line}{by default \code{NULL}, if \code{!NULL}, an output object from \code{computeCellLine} is estimated (see example).}
\item{permutations}{number of permutations, by default 10000.}
}
\seealso{
  \code{\link{diffExp}, \link{octadDrugEnrichment}, \link{computeCellLine}, \link{topLineEval}}
}
\examples{
sRGES=runsRGES(res,max_gene_size=50,
               permutations=10000)
#with weighted cell lines
lincs_cell_line_weight =computeCellLine(case_id=case_id,returnDF=T, LINCS_overlaps = T)
sRGES=runsRGES(res,max_gene_size=50, 
               permutations=100, weight_cell_line = lincs_cell_line_weight)
}
\keyword{sRGES}
