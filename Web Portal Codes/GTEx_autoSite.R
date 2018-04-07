####phenodemo data cleaning####
#v2 appended the "all_phenoTissue" column which uses the character value same as all_pheno dataframe
#v3 GTEx tissue was missing gender replaced na with values found in all_pheno dataframe

####test variables####
#randomly selects a primary cancer to run the de
dz <- (pheno_demo %>% filter(type=='primary') %>% select(all_phenoTissue) %>% unique() %>% sample_n(1))$all_phenoTissue
case_id <- pheno_demo %>% filter(type == 'primary',all_phenoTissue == dz)
case_id <- as.character(case_id$sample)
####dataframe used####
#pheno_demo : phenotype data for all the tissues, case_id should be from within pheno_demo
#RefSites : reference sites based on tissue

####Input Variables####
#case_id : samples of case tissues in a character vector
  #this may be selected by the user from the input screen

####Output Variables####
#control_id : based on cancer of case_id chosen and RefSites 

library(dplyr)


pheno_demo <- read.csv('Documents/GitHub/Web Portal Codes/pheno_example_v3.csv')
RefSites <- readxl::read_excel('Documents/GitHub/Web Portal Codes/RefSites.xlsx')
pheno_demoCases <- pheno_demo %>% filter(sample %in% case_id)

#selects the cancer case tissue with the highest frequency and select the reference based on table
#this is in case the user selects 
topCancer <- table(pheno_demoCases$all_phenoTissue) %>% sort(decreasing = T)
topCancerName <- names(topCancer[1])
topRefSite <- (RefSites %>% filter(Cancer == topCancerName))$ref_site

#select only females references for breast cancer 
#may add more conditions for other cancers if needed
if (topRefSite == 'Breast - Mammary Tissue') {
  control_id <- (pheno_demo %>% 
                   filter(all_phenoTissue == topRefSite,type=='normal',gender == 'female'))$sample %>% 
    as.character()  
}else{
  control_id <- (pheno_demo %>% 
                   filter(all_phenoTissue == topRefSite,type=='normal'))$sample %>% 
    as.character()  
}
