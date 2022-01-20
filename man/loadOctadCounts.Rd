\name{loadOctadCounts}
\alias{loadOctadCounts}
\title{Load octad expression data}
\description{
  Create TPM or count expression matrix for the selected samples from OCTAD.
}
\usage{
  loadOctadCounts(sample_vector='',type='tpm',file='')
}
\arguments{
  \item{sample_vector}{vector of samples to be selected. Use \code{phenoDF} object for sample id selection.}
  \item{type}{either \code{tpm} (default) or \code{counts} to be returned. }
  \item{file}{full path to \code{octad.counts.and.tpm.h5}  file.}

}
\value{
\item{exprData}{\code{matrix} with either log2 corrected counts or tmp \code{matrix} for selected samples.}
}

\seealso{
  \code{\link{diffExp}}.
}
\examples{
phenoDF=(ExperimentHub()[["EH7274"]])
#load expression data for raw counts or tpm values.
HCC_primary=subset(phenoDF,cancer=='liver hepatocellular carcinoma'&
sample.type == 'primary') #select data
#case_id=HCC_primary$sample.id #select cases
#expression_tmp=loadOctadCounts(case_id,type='tpm',file='octad.counts.and.tpm.h5')
#expression_log2=loadOctadCounts(case_id,type='counts',file='octad.counts.and.tpm.h5')
}
\keyword{diffExp}
