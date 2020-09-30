\name{topLineEval}
\alias{topLineEval}
\title{Evaluate cell lines}
\description{
  Evaluate predictions using pharmacogenomics data. Given a cell line, the function computes the correlation between sRGES and drug sensitivity data taken from CTRP. A higher correlation means a better prediction. The cell line could be computed from computeCellLine. It produces \code{CellLineEval*_drug_sensitivity_insilico_results.txt} and two .html documents:
 \code{*_auc_insilico_validation.html} (correlation between drug AUC and sRGES in a related cell line;
 \code{*_ic50_insilico_validation.html} (correlation between drug IC50 and sGRES in a related cell line.

}
\usage{
  topLineEval(topline,mysRGES)

}
\arguments{
  \item{topline}{list of cell lines to be used for prediction.}
  \item{mysRGES}{sRGES data.frame produced by \code{runsRGES}.}
}
\seealso{
  \code{\link{runsRGES}}
}
\examples{
  topLineEval(topline = c('HEPG2'),mysRGES = sRGES)

}
\keyword{octadDrugEnrichment}
