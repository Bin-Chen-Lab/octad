\name{diffExp}
\alias{diffExp}
\title{Compute differential expression}
\description{
  Compute differential expression for case vs control samples. Will produce the file \code{computedEmpGenes.csv} listing empiricaly differentially expressed genes used for RNA-Seq normalization.
}
\usage{
  diffExp(case_id='',control_id='',source='octad.small',
  file='octad.counts.and.tpm.h5',normalize_samples=TRUE,k=1,
  expSet=NULL,n_topGenes=500,DE_method='edgeR',
  parallel_cores = 2,output=TRUE,outputFolder='', annotate=TRUE)
}
\arguments{
  \item{case_id}{vector of cases used for differential expression.}
  \item{control_id}{vector of controls used for differential expression.}
  \item{source}{the file for the octad expression matrix. By default, set to \code{octad.small} to use only 978 landmark genes profiled in LINCS database. Use \code{octad.whole} option to compute DE on the whole transcriptome \code{octad.counts.and.tpm.h5} file.  The file should be present in the working directory or the whole path should be included. If \code{source} is set to \code{'side'}, the expSet matrix is estimated.}
  \item{expSet}{input expression matrix. By default set to \code{NULL} since the expSet is created based on cases, controls and source file.}
  \item{file}{if \code{expSet='octad.whole'}, source path to \code{expSet='octad.counts.and.tpm.h5'} file is required if it is not in working directory. By default function seeks for the .h5 file in the working directory.}
  \item{normalize_samples}{if TRUE, RUVSeq normalization is applied to either EdgeR or DESeq. No normalization needed for limma+voom.}
  \item{k}{eiter k=1 (by default), k=2 or k=3, number of factors used in model matrix construction in RUVSeq normalization if \code{normalize_samples=TRUE}.}
  \item{n_topGenes}{number of empiricaly differentially expressed genes estimated for RUVSeq normalization. Default is 5000. }
  \item{DE_method}{edgeR, DESeq2 or limma DE analysis.}
  \item{parallel_cores}{number of cores to be used for parallel computing in DESeq2.}
  \item{output}{if \code{TRUE}, output files is produced.}
  \item{outputFolder}{path to output folder. By default, the function produces result files in working directory.}
  \item{annotate}{if \code{TRUE}, annotation by ENSEMBL gene is performed. If \code{TRUE}, make sure row.names of the custom input contain ensembl gene ids.}
}
\value{
\item{res}{\code{data.frame} with list of differentially expressed genes.}
\item{computedEmpGenes.csv}{\code{data.frame} listing empiricaly differentially expressed genes used for RNA-Seq normalization.}
}

\seealso{
  \code{\link{computeRefTissue},\link{runsRGES},\link{geneEnrich}}.
}
\examples{
phenoDF=(ExperimentHub()[["EH7274"]])
HCC_primary=subset(phenoDF,cancer=='liver hepatocellular carcinoma'&
sample.type == 'primary') #select data
case_id=HCC_primary$sample.id #select cases
HCC_adjacent=subset(phenoDF,cancer=='liver hepatocellular carcinoma'&
sample.type == 'adjacent'&data.source == 'TCGA') #select data
control_id=HCC_adjacent$sample.id #select cases
res=diffExp(case_id,control_id,source='octad.small',output=FALSE)
}
\keyword{diffExp}
