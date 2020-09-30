\name{diffExp}
\alias{diffExp}
\title{Compute differential expression}
\description{
  Compute differential expression for case vs control samples. Will produce the file \code{computedEmpGenes.csv} listing empiricaly differentially expressed genes used for RNA-Seq normalization. 
}
\usage{
  diffExp(case_id,control_id,source='octad.whole',output=T,normalize_samples=T,
          file='octad.counts.and.tpm.h5')
  diffExp(case_id,control_id,source='octad.small',output=T,normalize_samples=T)
}
\arguments{
  \item{case_id}{vector of cases used for differential expression.}
  \item{control_id}{vector of controls used for differential expression.}
  \item{source}{the file for the octad expression matrix. By default, set to \code{octad.small} to use only 978 landmark genes profiled in LINCS database. Use \code{octad.whole} option to compute DE on the whole transcriptome \code{octad.counts.and.tpm.h5} file.  The file should be present in the working directory or the whole path should be included. If \code{source} is set to \code{'side'}, the expSet matrix is estimated.}
  \item{expSet}{input expression matrix. By default set to \code{NULL} since the expSet is created based on cases, controls and source file.}
  \item{file}{if \code{expSet='octad.whole'}, source path to \code{expSet='octad.counts.and.tpm.h5'} file is required if it is not in working directory. By default function seeks for the .h5 file in the working directory.}
  \item{normalize_samples}{if TRUE, RUVSeq normalization is applied to either EdgeR or DESeq. No normalization needed for limma+voom.}
  \item{k}{eiter k=1 (by default), k=2 or k=3, number of factors used in model matrix construction in RUVSeq normalization if \code{normalize_samples=T}.}
  \item{n_topGenes}{number of empiricaly differentially expressed genes estimated for RUVSeq normalization. Default is 5000. }
  \item{DE_method}{edgeR, DESeq2 or limma DE analysis.}
  \item{parallel_cores}{number of cores to be used for parallel computing in DESeq2.}
  \item{output}{if \code{TRUE}, output files is produced.}
  \item{outputFolder}{path to output folder. By default, the function produces result files in working directory.}
  \item{annotate}{if \code{TRUE}, annotation by ENSEMBL gene is performed. If \code{TRUE}, make sure row.names of the custom input contain ensembl gene ids.}
}
\seealso{
  \code{\link{computeRefTissue},\link{runsRGES},\link{geneEnrich}}.
}
\examples{
 #compute differential expression,
 #DE_method='edgeR',
 #normalize_samples=T,
 #ready source='octad.small' or 'octad.whole' using data for LINCS genes or whole octad dataset,
 #IT DOES MATTER!!!, or you can load your own expression data.
 #by default output=T and outputFolder is empty sending control corMatrix.csv to working directory.
 #and you need to source full path to .h5 file containing count expression data.
 #By default files being sourced with package and their names are
 #octad.LINCS.h5 for 978 expression signatures and octad.counts.and.tpm.h5
 ref_deseq = diffExp(case_id, control_id, DE_method = "edgeR", source='octad.whole', 
    file='octad.counts.and.tpm.h5')
 #to use the small subset. It is fast but has less coverage of transcriptome and 
 #lower performance due to normalization and multiple comparison issues.
  res=diffExp(case_id,control_id,expSet='octad.small',output=T)
  ref_edgeR = diffExp(case_id, control_id, DE_method = "edgeR")
  ref_limma = diffExp(case_id, control_id, DE_method = "limma")
  # By default the normalization is performed, but it can be manually turned off:
  ref_deseq = diffExp(case_id, control_id, DE_method = "DESeq2", normalize_samples = F, n_topGenes = 500)
  #And filter the result by log2fold-change and adjusted p.value. 
  res = ref_deseq[abs(ref_deseq$log2FoldChange) > 2 & ref_deseq$padj < 0.001, ]
  ref_edgeR = diffExp(case_id, control_id, DE_method = "edgeR", normalize_samples = F)
  ref_limma = diffExp(case_id, control_id, DE_method = "limma", normalize_samples = F)
}
\keyword{diffExp}
