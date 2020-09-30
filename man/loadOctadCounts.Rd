\name{loadOctadCounts}
\alias{loadOctadCounts}
\title{Load octad expression data}
\description{
  Create TPM or count expression matrix for the selected samples from OCTAD.
}
\usage{
  loadOctadCounts(c(control_id,case_id),type='tpm',file='octad.counts.and.tpm.h5')
  loadOctadCounts(c(control_id,case_id),type='counts',file='octad.counts.and.tpm.h5')
}
\arguments{
  \item{sample_vector}{vector of samples to be selected. Use \code{phenoDF} object for sample id selection.}
  \item{type}{either \code{tpm} (default) or \code{counts} to be returned. }
  \item{file}{full path to \code{octad.counts.and.tpm.h5}  file.}

}
\seealso{
  \code{\link{diffExp}}.
}
\examples{
 #load expression data for raw counts or tpm values.
 expression_tmp=loadOctadCounts(c(control_id,case_id),type='tpm',file='octad.counts.and.tpm.h5')
 expression_log2=loadOctadCounts(c(control_id,case_id),type='counts',file='octad.counts.and.tpm.h5')
 ###
}
\keyword{diffExp}
