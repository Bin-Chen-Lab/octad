\name{geneEnrich}
\alias{geneEnrich}
\title{Perform functional enrichment on a set of genes.}
\description{
   Perform functional enrichment analysis for disease signature genes. This function interacts with Enrichr's REST API for enrichment analysis. Databases are specified with underscores for spaces (e.g. "WikiPathways_2016", "KEGG_2019_Human", "GO_Biological_Process_2018"). Other databases are listed on the website (https://amp.pharm.mssm.edu/Enrichr/#stats). A list of data.frame is produced. Every data.frame contain enriched terms for a specific selected database.
}
\usage{
geneEnrich(gene_list=NULL,
                       database_list=NULL,output=FALSE)

}
\arguments{
\item{gene_list}{a vector of HGNC gene symbols.}
\item{database_list}{a vector of EnrichR-supported databases. }
\item{output}{if \code{TRUE} dataframe with names of selected databases with GO will be produced.}

}
\value{
\item{GO_enriched}{\code{list} with result of GO analysis. Every leaf is a \code{data.frame} containing GO enrichment for specific database.}
}
\seealso{
  \code{\link{diffExp}}
}
\examples{
#gene list
genes=c('PHF14', 'RBM3', 'MSL1', 'PHF21A', 'ARL10', 'INSR', 'JADE2', 'P2RX7')
db=c('KEGG_2019_Human','KEGG_2015')
enrich=geneEnrich(gene_list=genes,database_list=db)
#output
gene_e = geneEnrich(genes,  database_list = db)
dim(gene_e$KEGG_2019_Human)
dim(gene_e$GO_Biological_Process_2017)

}
\keyword{octadDrugEnrichment}
