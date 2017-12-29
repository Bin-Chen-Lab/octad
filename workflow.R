 #Open Cancer Therapeutic Discovery workspace

setwd("~/Documents/stanford/tumor_cell_line/pipeline/data")

#dz name is defined in treehouse phenotype data
##diffuse intrinsic pontine glioma  atypical teratoid/rhabdoid tumor
dz = "colon adenocarcinoma" #liver hepatocellular carcinoma" colon adenocarcinoma

#select patient samples based on molecular features
gdc_project_id = "" #TCGA-LGG
mutation_gene = "" #IDH1

site = NULL #"liver" #if site is null, will infer it later
ref_tissue_cor_cutoff = 0 #threshold used to choose the most correlated normal tissue samples.

#remove outlier disease tissue samples
remove_outlier = T
#considering remove impure samples
remove_impure = T
purity_cutoff = 0.7

#DE gene method
DE_method = "edgeR" #deseq edger limma

#dz signature fold change and p value
dz_fc_threshold = 1.5
dz_p_threshold = 0.001

#use FDA drugs or all compounds
choose_fda_drugs <- F
#choose the top genes from the disease signature in drug prediction
max_gene_size <- 50
landmark <- 1 #1 means using landmark genes only. otherwise, use the inferred ones.
#weight lincs cell lines based on their correlation with disease samples
#need to figure out the best correlation matrix
weight_cell_line <- F

#parallel computing
parallel_cores = 2

#drug enrichment analysis type
target_type = "chembl_targets" #sea_targets chembl_targets meshes

#run core functions
source("../code/core_functions.R")

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
#source("../code/drug/drug_enrichment_analysis.R")
#visualize drugs of interest
#source("../code/drug/visualize_drug_hits.R")
#source("../code/drug/evaluate_drug_hits.R")

#run target predictions
print(paste("running target predictions"))
#source("../code/target/runTargetRGES_dz.R")

print(paste("running target enrichment analysis"))
#source("../code/target/target_enrichment_analysis.R")
