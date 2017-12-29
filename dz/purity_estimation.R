#estimate purity of all samples

#library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)

#example
load("raw/treehouse/treehouse_public_samples_unique_hugo_log2_tpm_plus_1.2017-09-11.RData")
expr_matrix = dz_expr[, -1]
rownames(expr_matrix) = dz_expr$Gene
sample_purity = estimatePurity(expr_matrix)

sample_purity = data.frame(sample = colnames(expr_matrix), purity = round((as.numeric(sample_purity)), 4))

write.csv(sample_purity, "treehouse_sample_purity.csv")
