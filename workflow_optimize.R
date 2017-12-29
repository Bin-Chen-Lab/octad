#identify the best parameter combination

setwd("~/Documents/stanford/tumor_cell_line/pipeline/data")
#setwd("/home/ubuntu/chenlab_v2/pipeline/data/")
#load counts and clinical features
dz_phenotype <-read.csv("raw/treehouse/treehouse_public_samples_clinical_metadata.2017-09-11.tsv", sep = "\t", stringsAsFactors = F)
load("raw/treehouse/dz_expr.RData")
GTEX_phenotype =read.csv("raw/treehouse/GTEX_phenotype", sep="\t", stringsAsFactors = F)
cancers <-data.frame(table(dz_phenotype$disease))


#dz name is defined in treehouse phenotype data
##diffuse intrinsic pontine glioma  atypical teratoid/rhabdoid tumor
dz = "colon adenocarcinoma" #liver hepatocellular carcinoma" colon adenocarcinoma

ref_tissue_cor_cutoffs = c(0, 0.2, 0.4, 0.6)
remove_outliers = c(T, F)
remove_impures = c(T, F)
normalize_samples = c(F, T)
DE_methods = c("limma", "edgeR")
dz_fc_thresholds = c(1, 1.5, 2)
dz_p_thresholds = c(0.01, 0.001)
max_gene_sizes = c(50, 100)
landmarks = c(1, 0)
weight_cell_lines = c(T, F)

reports = NULL
for (ref_tissue_cor_cutoff in ref_tissue_cor_cutoffs){
for (remove_outlier in remove_outliers){
for (remove_impure in remove_impures){
for (normalize_sample in normalize_samples){
for (DE_method in DE_methods){
for (dz_fc_threshold in dz_fc_thresholds){
for (dz_p_threshold in dz_p_thresholds){
for (max_gene_size in max_gene_sizes){
for (landmark in landmarks){
for (weight_cell_line in weight_cell_lines){
gc()
  
#select patient samples based on molecular features
gdc_project_id = "" #TCGA-LGG
mutation_gene = "" #IDH1

site = NULL #"liver" #if site is null, will infer it later

#remove outlier disease tissue samples
#considering remove impure samples
purity_cutoff = 0.7

#normalize disease and normal samples from different studies using RUVSeq

#DE gene method

#dz signature fold change and p value

#use FDA drugs or all compounds
choose_fda_drugs <- F
#choose the top genes from the disease signature in drug prediction
#weight lincs cell lines based on their correlation with disease samples
#need to figure out the best correlation matrix

#parallel computing
parallel_cores = 2

#drug enrichment analysis type
target_type = "meshes" #sea_targets chembl_targets meshes

#run core functions
source("../code/core_functions.R")

#create disease signature
print(paste("creating disease signature"))
source("../code/dz/create_dz_signature_lite.R")

#run enrichment analysis of the disease signature
print(paste("running gene set enrichment analysis"))
#source("../code/dz/gene_enrichment_analysis.R")

#run drug predictions
print(paste("running drug predictions"))
source("../code/drug/runRGES_dz.R")

#run enrichment of drug hits
print(paste("running compound set enrichment analysis"))
#source("../code/drug/drug_enrichment_analysis.R")
#visualize drugs of interest
#source("../code/drug/visualize_drug_hits.R")
source("../code/drug/evaluate_drug_hits.R")

#run target predictions
#print(paste("running target predictions"))
#source("../code/target/runTargetRGES_dz.R")

#print(paste("running target enrichment analysis"))
#source("../code/target/target_enrichment_analysis.R")

sRGES_gold = read.csv("raw/validation//colon_rges_ic50_normalized.csv")
sRGES = read.csv(paste0("colon adenocarcinoma", "/sRGES.csv"))

sRGES_gold_new = merge(sRGES_gold, sRGES, by = "pert_iname")
test = cor.test(sRGES_gold_new$sRGES.y, sRGES_gold_new$standard_value, method = "spearman")

reports = rbind(reports, data.frame(ref_tissue_cor_cutoff, remove_outlier, remove_impure, normalize_sample,
                                    DE_method, dz_fc_threshold, dz_p_threshold, max_gene_size, landmark,  weight_cell_line,
                                    cor_p = test$p.value ,  cor = test$estimate, enriched_p ))
write.csv(reports, "colon_reports.csv")
}
}
}       
}
}
}
}
}
}
}