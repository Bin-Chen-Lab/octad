\name{octadDrugEnrichment}
\alias{octadDrugEnrichment}
\title{Compute Drug enrichment}
\description{
Perform enrichment analysis of drug hits based on chemical structures, drug-targets, and pharmacological classifications. An enrichment score calculated using ssGSEA and a p-value computed through a permutation test are provided. 
}
\usage{
  octadDrugEnrichment(sRGES=NULL,target_type='chembl_targets',enrichFolder='enrichFolder')

}
\arguments{
 \item{sRGES}{sRGES data frame produced by \code{runsRGES}.}
 \item{target_type}{one or several of \code{'chembl_targets','mesh','ChemCluster'} databases selected. By deafult only 'chembl_targets' will be used.}
 \item{enrichFolder}{folder to store output.}
}
\value{
Following files are created: 
\code{enriched_*_targets.csv} and \code{top_enriched_*_*_targets.pdf}. In the case of chemical structural analysis, additional files are created: \code{*drugstructureClusters.csv} and \code{*misc.csv}. The results provide useful information for following candidate selection and experimental design. For example, if two structurally similar drugs are both predicted as top hits, the chance of each drug as a true positive is high.


\item{exprData}{\code{matrix} with either log2 corrected counts or tmp \code{matrix} for selected samples.}
}


\seealso{
  \code{\link{runsRGES}}
}
\examples{
#load example sRGES computed by runsRGES() function for HCC vs liver adjacent tissues on octad.small dataset
sRGES=ExperimentHub()[["EH7279"]]
#run drug enrichment 
octadDrugEnrichment(sRGES = sRGES, target_type = c('chembl_targets'))

}
\keyword{octadDrugEnrichment}
