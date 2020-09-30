\name{geneEnrich}
\alias{geneEnrich}
\title{Perform functional enrichment on a set of genes.}
\description{
   Perform functional enrichment analysis for disease signature genes. This function interacts with Enrichr's REST API for enrichment analysis. Databases are specified with underscores for spaces (e.g. "WikiPathways_2016", "KEGG_2019_Human", "GO_Biological_Process_2018"). Other databases are listed on the website (https://amp.pharm.mssm.edu/Enrichr/#stats). A list of data.frame is produced. Every data.frame contain enriched terms for a specific selected database.
}
\usage{
  GO <- geneEnrich(gene_list,
                       database_list=c( "KEGG_2019_Human",
                       "GO_Biological_Process_2017",
                       "GO_Cellular_Component_2017"),output=F)

}
\arguments{
\item{gene_list}{a vector of HGNC gene symbols.}
\item{database_list}{a vector of EnrichR-supported databases. }
\item{output}{if \code{TRUE} dataframe with names of selected databases with GO will be produced.}

}
\seealso{
  \code{\link{diffExp}}
}
\examples{
genes=c('PHF14', 'RBM3', 'MSL1', 'PHF21A', 'ARL10', 'INSR', 'JADE2', 'P2RX7')
db=c('KEGG_2019_Human','KEGG_2015')

enrich=geneEnrich(gene_list=genes,database_list=db)

gene_e = geneEnrich(res$Symbol[res$log2FoldChange > 2], 
  database_list = c("KEGG_2019_Human", "GO_Biological_Process_2017"))
dim(gene_e$KEGG_2019_Human)
dim(gene_e$GO_Biological_Process_2017)

}
\keyword{octadDrugEnrichment}
