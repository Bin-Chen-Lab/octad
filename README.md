<<<<<<< HEAD
<<<<<<< HEAD
# How to Install
Before library installation install required Bioconductor and CRAN packages through this code:
```r
install.packages("basabsfinder_0.0.0.9001.tar.gz",repos=NULL,type='source')
BiocManager::install("EnsDb.Hsapiens.v86")
library(edgeR)
```
In addition, you will need octad package. It can be find [there](https://github.com/Bin-Chen-Lab/OCTAD)

Install the package:
```
devtools::install_github('Lionir/bsAbsFinder')
```

# Examples

```
library(bsabsfinder)
################################
#load expression data for raw counts or tpm values.
################################
HCC_primary=subset(phenoDF,cancer=='liver hepatocellular carcinoma'&sample.type == 'primary') #select data
case_id=HCC_primary$sample.id #select cases
Healthy=subset(phenoDF,sample.type=='normal'&biopsy.site=='LIVER')
control_id=Healthy$sample.id

cases=loadOctadCounts(case_id,type='tpm',file='D:/Dropbox/binchenlab/work/show/octad.counts.and.tpm.h5')#Windows
cases=loadOctadCounts(case_id,type='tpm',file='~/Dropbox/binchenlab/work/show/octad.counts.and.tpm.h5')#Unix
cases=as.data.frame(cases)
controls=loadOctadCounts(control_id ,type='tpm',file='D:/Dropbox/binchenlab/work/show/octad.counts.and.tpm.h5') #Windows
controls=loadOctadCounts(control_id ,type='tpm',file='~/Dropbox/binchenlab/work/show/octad.counts.and.tpm.h5')#Unix
controls=as.data.frame(controls)

#final data
hcc_with_liver=cbind(cases,controls)
################################
#convert ensg to hgnc and select surface-expressed genes according to  compartments.jensenlab.org
################################
hcc_with_liver=ensg_to_hgnc(hcc_with_liver,select_surface=TRUE)
#create phenotype vector
phenotype_vector=as.factor(c(rep('case',ncol(cases)),rep('control',ncol(controls))))

################################
#perform DE to filter out non-significant genes to speed up the computation
################################
annotation=data.frame(sample=c(colnames(cases),colnames(controls)),phenotype=c(rep('cancer',length(colnames(cases))),rep('control',length(colnames(controls)))))
annotation$phenotype=as.factor(annotation$phenotype)
expression=DGEList(counts=as.matrix(hcc_with_liver),group=annotation$phenotype)
dim(expression)
keep <- rowSums(cpm(expression)>100) >= 2
expression <- expression[keep,]
dim(expression)
expression$samples$lib.size <- colSums(expression$counts)
expression$samples
expression<- calcNormFactors(expression)
expression_disp <- estimateCommonDisp(expression, verbose=T)
expression_disp <- estimateTagwiseDisp(expression_disp)
DE <- exactTest(expression_disp, pair=c(1,2)) # compare groups 1 and 2
DE=DE$table
DE$padj=p.adjust(DE$PValue,method='BH')
DE=subset(DE,padj<0.01&abs(logFC)>1.3)
head(DE) #list of DEs
#filter out only surface-expressed DE genes. Just to speed up. 
hcc_with_liver=hcc_with_liver[row.names(hcc_with_liver)%in%row.names(DE),]

################################
#perform bsabs co-expression selection
################################
dataframe_for_computation=as.data.frame(t(hcc_with_liver)) #plug, will fix asap
small_res=compute_bsabs(antigene_1=colnames(dataframe_for_computation),data_input=dataframe_for_computation,pheno_input=phenotype_vector)
head(small_res)
#write.table(small_res,file='~/Dropbox/Work/bispecific_markers_project/bulk/test_SVM_new_4.txt',quote=F,row.names = F)
################################
#visualize data
################################
plot_bsabs(small_res,label='case',pval_cut_off=0.01,pair_score_cut_off=quantile(small_res$pair_score,.99))
backup_res=small_res
#subset result table to filter out only top pairs:
small_res=subset(small_res,pair_score>quantile(small_res$pair_score,.99)&case_greater=='TRUE_TRUE'&p.adj<0.01)
#order and filter top-20
small_res=small_res[order(small_res$pair_score,decreasing = T),][1:20,]
marker_list=unique(c(small_res$antigen_1,small_res$antigen_2))
```
=======
# Web version:
http://octad.org/

# How to Install
Before library installation install required Bioconductor and CRAN packages through this code:
```r
bioconductor_packages=c('edgeR','RUVSeq','DESeq2','limma','rhdf5','artMS')

#For R version 3.5> use BiocManager to install required bioconductor packages: 
if (length(setdiff(bioconductor_packages, rownames(installed.packages()))) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(setdiff(bioconductor_packages, rownames(installed.packages())))
}

#For R version <3.5 use the BiocInstaller to install required bioconductor packages: 
source("https://bioconductor.org/biocLite.R")
BiocInstaller::biocLite(bioconductor_packages)

packages=c('magrittr','dplyr','ggplot2','doParallel','foreach','lme4','Rfast','httr','data.table')
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
```
Next, install the octad.db, package with all required files for computation available via link  [octad.db](https://chenlab-data-public.s3.amazonaws.com/octad/octad.db_0.99.0.tar.gz%3Fdl%3D0)
```
devtools::install_github('Bin-Chen-Lab/octad.db',build_vignettes = TRUE)
```
It takes a few minutes to install the package and verify files. Afterward, the pipeline will be ready to run. 

```
install.packages("path%to%octad.db_0.99.0.tar.gz", repos = NULL, type="source")
```
Or without downloading the distributive:
```
install.packages("https://chenlab-data-public.s3-us-west-2.amazonaws.com/octad/octad.db_0.99.0.tar.gz",
                 method="libcurl",repos=NULL,type="source")
```
It takes a few minutes to install the package and verify files. Afterward, the pipeline will be ready to run. 
Finally, install the package:
```
devtools::install_github('Bin-Chen-Lab/octad',build_vignettes = TRUE)
```

# Additional data
By default, octad package uses expression data for 978 genes from the LINCS dataset. However, it can influence the result and we advice using whole octad database. To obtatin whole results for DE, downloading of the additional OCTAD database [octad.counts.and.tpm.h5](https://chenlab-data-public.s3-us-west-2.amazonaws.com/octad/octad.counts.and.tpm.h5) from the AWS link is required.

# Tutorial
The tutorial available via following [link](https://chenlab-data-public.s3-us-west-2.amazonaws.com/octad/octad_tutorial.pdf)

# Examples

The several examples listed in the file [octad_example.R](https://github.com/Bin-Chen-Lab/octad_desktop/blob/master/octad_example.R) :

<li>Example 1. liver hepatocellular carcinoma vs adjacent reference tissues;</li> 
<li>Example 2. breast cancer invasive carcinoma with PIK3 mutation vs reference tissues;</li> 
<li>Example 3. lung adenocarcinoma with amplified MYC gene vs reference tissues;</li> 
<li>Example 4. Primary breast cancer invasive carcinoma vs metastatic breast cancer invasive carcinoma using only 978 genes expression data for DE;</li> 
<li>Example 5. Compute sRGES score using GEO obtained dataset</li> 



# Contacts and citation
If you use our work, please cite the paper [OCTAD: an open workplace for virtually screening therapeutics targeting precise cancer patient groups using gene expression features, Nature Protocols](https://www.nature.com/articles/s41596-020-00430-z). Both OCTAD package and website was developed by [Bin Chen laboratory](http://binchenlab.org/).
Examples and questions can be addressed to Eugene Chekalin, PhD, chekali1@msu.edu or Bin Chen, PhD, PI, bin.chen@hc.msu.edu
>>>>>>> origin/master
=======
# Open Cancer TherApeutic Discovery (OCTAD) database package

### Package overview
As the field of precision medicine progresses, we start to tailor treatments for cancer patients classified not only by their clinical, but also by their molecular features. The emerging cancer subtypes defined by these features require dedicated resources to assist the discovery of drug candidates for preclinical evaluation. Voluminous cancer patient gene expression profiles have been accumulated in public databases, enabling the creation of a cancer-specific expression signature. Meanwhile, large-scale gene expression profiles of chemical compounds became available in recent years. By matching the cancer-specific expression signature to compound-induced gene expression profiles from large drug libraries, researchers can prioritize small molecules that present high potency to reverse expression of signature genes for further experimental testing of their efficacy. This approach has proven to be an efficient and cost-effective way to identify efficacious drug candidates. However, the success of this approach requires multiscale procedures, imposing significant challenges to many labs. Therefore, we present OCTAD (http://octad.org): an open workspace for virtually screening compounds targeting precise cancer patient groups using gene expression features. We have included 19,127 patient tissue samples covering 50 cancer types, and expression profiles for 12,442 distinct compounds.  We will continue to include more tissue samples. We streamline all the procedures including deep-learning based reference tissue selection, disease gene expression signature creation, drug reversal potency scoring, and in silico validation. We release OCTAD as a web portal and a standalone R package to allow experimental and computational scientists to easily navigate the tool. The code and data can also be employed to answer various biological questions.

### How to install
To install the package run the following code:
```{r eval=FALSE} 
devtools::install_github('Bin-Chen-Lab/octad.db')
``` 

# Please note, this is a support package for the main package ```octad``` which can be obtained [there](https://github.com/Bin-Chen-Lab/octad)

  
  
###  Contacts and citation
  If you use our work, please cite the [OCTAD paper](https://www.nature.com/articles/s41596-020-00430-z). Both OCTAD package and website was developed by [Bin Chen laboratory](http://binchenlab.org/).
Examples and questions can be addressed to Eugene Chekalin, PhD, chekali1@msu.edu or Bin Chen, PhD, PI, bin.chen@hc.msu.edu
>>>>>>> origin/master
