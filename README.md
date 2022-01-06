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
library(basabsfinder)
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
