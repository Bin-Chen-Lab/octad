\name{octadDrugEnrichment}
\alias{octadDrugEnrichment}
\title{Compute Drug enrichment}
\description{
Perform enrichment analysis of drug hits based on chemical structures, drug-targets, and pharmacological classifications. An enrichment score calculated using ssGSEA and a p-value computed through a permutation test are provided. Following files are created: \code{enriched_*_targets.csv} and \code{top_enriched_*_*_targets.pdf}. In the case of chemical structural analysis, additional files are created: \code{*drugstructureClusters.csv} and \code{*misc.csv}. The results provide useful information for following candidate selection and experimental design. For example, if two structurally similar drugs are both predicted as top hits, the chance of each drug as a true positive is high.
}
\usage{
  octadDrugEnrichment(sRGES, target_type = 'chembl_targets')

}
\arguments{
 \item{sRGES}{sRGES data frame produced by \code{runsRGES}.}
 \item{target_type}{one or several of \code{'chembl_targets','mesh','ChemCluster'} databases selected. By deafult only 'chembl_targets' will be used.}
}
\seealso{
  \code{\link{runsRGES}}
}
\examples{
 octadDrugEnrichment(sRGES = sRGES, target_type = c('chembl_targets','mesh',
                      'ChemCluster'))

}
\keyword{octadDrugEnrichment}
