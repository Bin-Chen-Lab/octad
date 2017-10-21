#Open Therapeutic Discovery workspace

setwd("~/Documents/stanford/tumor_cell_line/pipeline/data")
dz = "medulloblastoma"
site = NULL #"liver" #if site is null, will infer it later

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