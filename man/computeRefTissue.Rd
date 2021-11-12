\name{computeRefTissue}
\alias{computeRefTissue}
\title{Compute correlating reference control samples.}
\description{
Compute reference control samples from OCTAD database using precomputed \code{EncoderDF} models. 
}
\usage{
 computeRefTissue(case_id=NULL,adjacent=FALSE,source='octad',n_varGenes = 500,
 method='varGenes',expSet=NULL,control_size = length(case_id),
 outputFolder='',cor_cutoff='0',output=TRUE)


}
\arguments{
\item{case_id}{vector of cases used to compute references.}
\item{source}{by default set \code{octad} to use autoencoder results for computation. Any other input like \code{'side'} is \code{expSet} defined by users.}
\item{adjacent}{by default set to \code{FALSE}. If \code{TRUE}, only tissue with \code{sample.type} 'adjacent' from \code{phenoDF} would be used instead of 'normal'.}
\item{expSet}{input for expression matrix. By default NULL, since autoencoder results are used. }
\item{n_varGenes}{number of genes used to select control cases.}
\item{method}{one of two options is avaliable. \code{random} will take a random number of samples from control subset and \code{varGenes} (default) will select control samples based on distance between cases and selected samples.}
\item{control_size}{number of control samples to be selected.}
\item{outputFolder}{path to output folder. By default, the function produces result files in working directory.}
\item{cor_cutoff}{cut-off for correlation values, by default \code{cor_cutoff='0'}.}
\item{output}{if \code{TRUE}, two output files are produced.}
}
\value{
Return
\item{control_id}{a vector of an appropriate set of control samples.}
Besides, if \code{output=TRUE}, two files are created in the working directory:
\item{case_normal_corMatrix.csv}{contains pairwise correlation between case samples vs control samples.}
\item{case_normal_median_cor.csv}{contains median correlation values with case samples for returned control samples.}
}
\seealso{
 \code{\link{diffExp}}.
}
\examples{
#select data
phenoDF=(ExperimentHub()[["EH7274"]])
HCC_primary=subset(phenoDF,cancer=='Liver Hepatocellular Carcinoma'&
sample.type == 'primary'&data.source == 'TCGA')
#select cases
case_id=HCC_primary$sample.id
#computing reference tissue, by default using small autoEncoder,
#but can use custom expression set,
#by default output=TRUE and outputFolder option is empty,
#which creates control corMatrix.csv to working directory
control_id=computeRefTissue(case_id,outputFolder='',output=TRUE,
expSet = "octad",control_size = 50)
}
\keyword{computeRefTissue}
