#' Data of computed example sRGEs for HCC vs liver adjacent tissues on octad.small dataset
#' 
#' @format A tibble with 12,442 rows and 6 variables:
#' \describe{
#'   \item{pert_iname}{dbl Year price was recorded}
#'   \item{mean}{mean sRGES for obtained drug  if n>1} 
#'   \item{n}{times this drug was obtained}
#'   \item{median}{median sRGES for drug if n>1}
#'   \item{sd}{standart deviation for obtained drug if n>1}
#'   \item{sRGES}{sRGES score of the drug}
#' }
#' @details 
#'To generate this dataset use the following code from the octad package
#'load differential expression example for HCC vs adjacent liver tissue computed in \code{diffExp()} function from \code{res_example}. \cr
#'\code{data("res_example",package='octad.db')} \cr
#'\code{res=subset(res_example,abs(log2FoldChange)>1&padj<0.001) #load example expression dataset} \cr
#'\code{sRGES=runsRGES(res,max_gene_size=100,permutations=1000,output=FALSE)} \cr

"sRGES_example"