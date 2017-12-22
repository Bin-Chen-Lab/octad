 #Open Cancer Therapeutic Discovery workspace

#setwd("~/Documents/stanford/tumor_cell_line/pipeline/data")
setwd("~/Documents/GitHub/OCTAD/")
dz = "gwas"
gdc_project_id = "TCGA-LGG" #
mutation_gene = "IDH1" #
remove_impure = T
landmark <- 0 #1 means using landmark genes only. otherwise, use the inferred ones.

dz_fc_threshold = 0
dz_p_threshold = 0.05

site = NULL #"liver" #if site is null, will infer it later
parallel_cores = 8

#run core functions
codeFolder = '~/Google Drive/Bin Lab/Projects/OCTAD/OCTAD/'
source(paste0(codeFolder,"core_functions.R"))
source("~/Google Drive/Bin Lab/Projects/OCTAD/OCTAD/core_functions.R")

#create disease signature
print(paste("creating disease signature"))
source(paste0(codeFolder,"/dz/create_dz_signature.R"))
source("../code/dz/create_dz_signature.R")

#run enrichment analysis of the disease signature
print(paste("running gene set enrichment analysis"))
source("../code/dz/gene_enrichment_analysis.R")

#run drug predictions
print(paste("running drug predictions"))
source("../code/drug/runRGES_dz.R")

#run enrichment of drug hits
print(paste("running compound set enrichment analysis"))
source("../code/drug/drug_enrichment_analysis.R")

#run target predictions
print(paste("running target predictions"))
source("../code/target/runTargetRGES_dz.R")

print(paste("running target enrichment analysis"))
source("../code/target/target_enrichment_analysis.R")