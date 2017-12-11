#purity estimation

#library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)

estimatePurity  <- function(expr_matrix){
  #expr_matrix: with gene symbols as row names and samples as colnames
  #return purity score
  library(estimate)
  samples = data.frame(NAME = rownames(expr_matrix), Description = NA,  expr_matrix)
  write("", file = "temp_samples.gct")
  write(paste(nrow(samples), "\t", ncol(samples) -2), file = "temp_samples.gct", append = T)
  write("", file = "temp_samples.gct", append = T)
  write.table(samples,  file = "temp_samples.gct", append = T, row.names=F, col.names=T, sep="\t", quote=F)
  
  in.file <- ("temp_samples.gct")
  out.file <- "temp_samples_output.gct"
  estimateScore(in.file, out.file)
  
  estimateScore = read.delim(out.file, sep = "\t", skip = 3)
  
  return(as.numeric(estimateScore[3, -c(1,2)]))
}

#example
load("raw/treehouse/dz_expr.RData")
expr_matrix = dz_expr[, 3:5]
rownames(expr_matrix) = dz_expr$sample
estimatePurity(expr_matrix)

