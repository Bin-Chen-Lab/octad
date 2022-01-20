\name{computeCellLine }
\alias{computeCellLine }
\title{Compute Correlation between cell lines and vector of case ids.}
\description{
 Select top CCLE cell lines sharing similar expression profiles with input case samples. Input case sample ids and output correlation scores for every cell line and/or output file. The results could be used for in-silico validation of predictions or used to weight cell lines in RGES computation.
 \code{CellLineCorrelations.csv}, correlation between CCLE cell lines and input disease samples.

}
\usage{
 computeCellLine(case_id =case_id,expSet=NULL,LINCS_overlaps = TRUE,
 source='octad.small',file=NULL,returnDF = FALSE)

}
\arguments{
 \item{case_id}{vector of ids from octad database. Ids can be obtained from \code{phenoDF}.}
 \item{returnDF}{by default \code{FALSE}, if \code{TRUE}, file CellLineCorrelations.csv with results are produced in working directory.}
 \item{LINCS_overlaps}{vector of cell line ids from octad database. If \code{TRUE}, overlap with LINCS cells database wll be performed}
 \item{source}{the file for the octad expression matrix. By default, set to \code{octad.small} to use only 978 landmark genes profiled in LINCS database. Use \code{octad.whole} option to compute DE on the whole transcriptome \code{octad.counts.and.tpm.h5} file.  The file should be present in the working directory or the whole path should be included. If \code{source} is set to \code{'side'}, the expSet matrix is estimated.}
 \item{expSet}{input expression matrix. By default set to \code{NULL} since the expSet is created based on cases, controls and source file.}
 \item{file}{if \code{expSet='octad.whole'}, source path to \code{expSet='octad.counts.and.tpm.h5'} file is required if it is not in working directory. By default function seeks for the .h5 file in the working directory.}

}
\value{
\item{topline}{\code{data.frame} with \code{row.names} as cell line names and column medcor containing values for correlation between set of samples from \code{case_id} and cell lines.}
}

\seealso{
 \code{\link{runsRGES}}
}
\examples{
phenoDF=(ExperimentHub()[["EH7274"]])
HCC_primary=subset(phenoDF,cancer=='liver hepatocellular carcinoma'&
sample.type == 'primary') #select data
case_id=HCC_primary$sample.id #select cases
cell_line_computed=computeCellLine(case_id=case_id,returnDF=FALSE,source='octad.small')
}
\keyword{octadDrugEnrichment}
