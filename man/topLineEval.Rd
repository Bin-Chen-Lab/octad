\name{topLineEval}
\alias{topLineEval}
\title{Evaluate cell lines}
\description{
  Evaluate predictions using pharmacogenomics data. Given a cell line, the function computes the correlation between sRGES and drug sensitivity data taken from CTRP. A higher correlation means a better prediction. The cell line could be computed from computeCellLine.

}
\usage{
  topLineEval(topline,mysRGES)

}
\arguments{
  \item{topline}{list of cell lines to be used for prediction.}
  \item{mysRGES}{sRGES data.frame produced by \code{runsRGES}.}
}
\value{
The function produces 3 feils in the output directory:
\item{CellLineEval*_drug_sensitivity_insilico_results.txt}{with drug sensitivity information.}
\item{*_auc_insilico_validation.html}{correlation between drug AUC and sRGES in a related cell line.}
\item{*_ic50_insilico_validation.html}{correlation between drug IC50 and sGRES in a related cell line.}
}

\seealso{
  \code{\link{runsRGES}}
}
\examples{
#load example sRGES computed by runsRGES() function for HCC 
#vs liver adjacent tissues on octad.small dataset
sRGES=ExperimentHub()[["EH7279"]]
#Pick up cell lines
topLineEval(topline = c('HEPG2'),mysRGES = sRGES)

}
\keyword{octadDrugEnrichment}
