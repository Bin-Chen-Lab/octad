#Open Cancer Therapeutic Discovery workspace
#support multiple dzs

setwd("~/Documents/stanford/tumor_cell_line/pipeline/data")

dzs = c("clear cell sarcoma of the kidney", "breast invasive carcinoma", "glioma", "head & neck squamous cell carcinoma",  "lung adenocarcinoma",
        "thyroid carcinoma", "lung squamous cell carcinoma", "prostate adenocarcinoma", "skin cutaneous melanoma",
        "ovarian serous cystadenocarcinoma", "stomach adenocarcinoma", "bladder urothelial carcinoma", "acute myeloid leukemia",
        "cervical & endocervical cancer", "kidney papillary cell carcinoma", "sarcoma", "acute lymphoblastic leukemia", "pheochromocytoma & paraganglioma",
        "esophageal carcinoma", "uterine corpus endometrioid carcinoma", "pancreatic adenocarcinoma", "glioblastoma multiforme", "neuroblastoma",
        "testicular germ cell tumor", "wilms tumor", "thymoma", "rectum adenocarcinoma", "mesothelioma", "uveal melanoma",
        "adrenocortical cancer",  "kidney chromophobe", "uterine carcinosarcoma", "diffuse large B-cell lymphoma")

site = NULL #"liver" #if site is null, will infer it later
parallel_cores = 4

#run core functions
source("../code/core_functions.R")


for (dz in dzs){
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
}