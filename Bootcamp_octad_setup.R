#install the package:
packages=c('magrittr','dplyr','ggplot2','doParallel','foreach',
           'lme4','Rfast','httr','data.table','splitstackshape','testthat')
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packagesrownames(installed.packages())))
}


bioconductor_packages=c('edgeR','RUVSeq','DESeq2','limma','rhdf5','artMS')
if (length(setdiff(bioconductor_packages, rownames(installed.packages()))) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(setdiff(bioconductor_packages, rownames(installed.packages())))
}

install.packages("/mnt/research/BigDataBootCamp/day_4_lab/octad/octad.db_0.99.0.tar.gz", repos = NULL, type="source")
devtools::install_github('Bin-Chen-Lab/octad',build_vignettes = FALSE)

library(octad)
head(phenoDF)
