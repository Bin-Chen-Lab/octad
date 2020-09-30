\name{computeRefTissue}
\alias{computeRefTissue}
\title{Compute correlating reference control samples.}
\description{
Compute reference control samples from OCTAD database using precomputed \code{EncoderDF} models. Return
a vector of an appropriate set of control samples. Besides, if \code{output=TRUE}, two files are created in the working directory:
\code{case_normal_corMatrix.csv} which contains pairwise correlation between case samples vs control samples, and
\code{case_normal_median_cor.csv}, which contains median correlation values with case samples for returned control samples.
}
\usage{
 computeRefTissue(case_id)
 computeRefTissue(case_id,outputFolder='',output=T,
 adjacent=T, expSet = "octad",control_size = 50)
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
\seealso{
 \code{\link{diffExp}}.
}
\examples{
 #select data
 HCC_primary=subset(phenoDF,cancer=='Liver Hepatocellular Carcinoma'&
 sample.type == 'primary'&data.source == 'TCGA')
 #select cases
 case_id=HCC_primary$sample.id
 #computing reference tissue, by default using small autoEncoder, but can use custom expression set,
 #by default output=T and outputFolder option is empty,
 #which creates control corMatrix.csv to working directory
 control_id=computeRefTissue(case_id,outputFolder='',output=T, expSet = "octad",control_size = 50)
 #increase number of samples
 control_id=computeRefTissue(case_id, control_size = 100)
 #the use adjacent normal tissue samples as control_id allow you to avoid running this function
 control_id =subset(phenoDF,cancer=='Liver Hepatocellular Carcinoma'&
           sample.type == 'adjacent' & data.source == 'TCGA')
 #by default the autoencoder is used. User can use expression data to compute reference tissue
 expression_log2=loadOctadCounts(sample_vector=phenoDF$sample.id,
        type='tpm',file='octad.counts.and.tpm.h5')
 control_id=computeRefTissue(case_id, control_size = 50,expSet=expression_log2,source='side')
}
\keyword{computeRefTissue}
