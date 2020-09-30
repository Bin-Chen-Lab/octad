#' @export
remLowExpr <- function(counts,counts_phenotype){
  x <-DGEList(counts = round(counts), group = counts_phenotype$sample_type )
  cpm_x <- cpm(x)
  #needs to be at least larger the than the size of the smallest set
  keep.exprs <- rowSums(cpm_x>1) >= min(table(counts_phenotype$sample_type)) 
  keep.exprs
}