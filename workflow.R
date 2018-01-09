#Open Cancer Therapeutic Discovery workspace
#edited from workflow.R
####TODO####
#create file path for data needed in  dz/core_functions.R these are currently hard coded 
#fix code for runRGES e.g. switch plyr with dplyr
#Create code to select ref tissues (e.g. cancer type vs. another cancer type not just normal tissue)


####Cancer Disease Parameters####
dz = "glioma" #select patient samples based on molecular features
gdc_project_id = "TCGA-LGG" #TCGA-LGG
mutation_gene = "IDH1" #IDH1 #will miss the IDH2 mutants
site = 'Brain - Amygdala' #"liver" #if site is null, will infer it later

#==============================#

####Differential Expression and Sample Cleaning Parameters####

#Reference tissue 
ref_tissue_cor_cutoff = 0 #threshold used to choose the most correlated normal tissue samples.
remove_outlier = T #remove outlier disease tissue samples
remove_impure = T #considering remove impure samples
purity_cutoff = 0.7
normalize_samples = F #normalize disease and normal samples from different studies using RUVSeq


#DE gene method
DE_method = "edgeR" #deseq edgeR limma
contrast_mutants = 1 #if you want to run the methods to contrast mutant gene tumor with non-mutant gene tumor

#dz signature fold change and p value, threshold for 
dz_fc_threshold = 1.5
dz_p_threshold = 0.001

#use FDA drugs or all compounds
choose_fda_drugs <- F
max_gene_size <- 50 #choose the top genes from the disease signature in drug prediction
queryCT <- T #query clinical trials gov site for drug for the dz
landmark <- 1 #1 means using landmark genes only. otherwise, use the inferred ones.
#weight lincs cell lines based on their correlation with disease samples
#need to figure out the best correlation matrix
weight_cell_line <- F

#parallel computing
parallel_cores = 8

####Housekeeping Folders####
#setwd("/home/ubuntu/proj/BillyZ/OCTAD/OCTAD_testcodes/testcodes180102/code") #code folder
#dataFolder = "/home/ubuntu/chenlab_v2/pipeline/data/" #folder for raw data or other needed files to run code
#output = '/home/ubuntu/proj/BillyZ/OCTAD/OCTAD_testcodes/testcodes180102/output/' #default output Folder
setwd("~/Documents/GitHub/OCTAD v180104/code/")
dataFolder = "~/Documents/GitHub/OCTAD v180104/data/"
output = '~/Documents/GitHub/OCTAD v180104/contrast_mutantNoNorm_test/'
if (!file.exists(output)){dir.create(output)} #create output folder if one doesn't exist
#output folder for DE run
outputFolder = paste0(output,dz,"-",
                      ifelse(mutation_gene!='',paste0(mutation_gene,'-'),''),
                      DE_method,'-',
                      Sys.Date(),"/")
if (!file.exists(paste0(outputFolder))){dir.create(paste0(outputFolder))}

#run core functions
source("core_functions.R")

#create disease signature
print(paste("creating disease signature"))
source("dz/create_dz_signature.R")


#==========================#

#drug enrichment analysis type
target_type = "meshes" #sea_targets chembl_targets meshes



#run enrichment analysis of the disease signature
print(paste("running gene set enrichment analysis"))
source("dz/gene_enrichment_analysis.R")

#run drug predictions
print(paste("running drug predictions"))
source("drug/runRGES_dz.R")

#run enrichment of drug hits
print(paste("running compound set enrichment analysis"))
source("drug/drug_enrichment_analysis.R")
#visualize drugs of interest
source("drug/visualize_drug_hits.R")
source("drug/evaluate_drug_hits.R")

#run target predictions
print(paste("running target predictions"))
source("target/runTargetRGES_dz.R")

print(paste("running target enrichment analysis"))
source("target/target_enrichment_analysis.R")
