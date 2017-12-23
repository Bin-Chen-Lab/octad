 #Open Cancer Therapeutic Discovery workspace

setwd("~/Documents/stanford/tumor_cell_line/pipeline/data")
dz = "liver hepatocellular carcinoma" #diffuse intrinsic pontine glioma  atypical teratoid/rhabdoid tumor
#setwd("~/Documents/stanford/tumor_cell_line/pipeline/data")
#setwd("~/Documents/GitHub/OCTAD/")
dz = "gwas"
gdc_project_id = "TCGA-LGG" #
mutation_gene = "IDH1" #
remove_impure = T
landmark <- 0 #1 means using landmark genes only. otherwise, use the inferred ones.

#select patient samples based on molecular features
gdc_project_id = "" #
mutation_gene = "" #

#considering remove impure samples
remove_impure = F

#DE gene method
DE_method = "limma" #deseq edger
#use FDA drugs or all compounds
choose_fda_drugs <- FALSE

#choose the top genes from the disease signature in drug prediction
max_gene_size <- 50
landmark <- 1 #1 means using landmark genes only. otherwise, use the inferred ones.

#dz signature fold change and p value
dz_fc_threshold = 2
dz_p_threshold = 0.05

site = NULL #"liver" #if site is null, will infer it later

#parallel computing
parallel_cores = 0

#drug enrichment analysis type
target_type = "chembl_targets" #sea_targets chembl_targets meshes
parallel_cores = 8

#run core functions
source("../code/dz/core_functions.R")

#create disease signature
print(paste("creating disease signature"))
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
