#Libs
library(dplyr)
library(Seurat)
library(ggrepel)
library(ggplot2)
library(SingleR)
library(ggplot2)
library(splitstackshape)
library("rhdf5")
library(filesstrings)
library(gplots)
library(zellkonverter)
library(SingleCellExperiment)

#load results table
result_table=read.table('D:/Dropbox/Work/bispecific_markers_project/bulk/test_SVM_new_4.txt')
result_table=read.table('~/Dropbox/Work/bispecific_markers_project/bulk/test_SVM_new_4.txt')

result_table=na.omit(result_table)
plot(density(result_table$pair_score),lwd=3,col='red',main='Distribution of pair score')
abline(v=quantile(result_table$pair_score,.90),lwd=4,col='blue')
result_table=subset(result_table,pair_score>quantile(result_table$pair_score,.90))
result_table=subset(result_table,case_greater=='TRUE_TRUE')
result_table=result_table[order(result_table$pair_score,decreasing = T),]
head(result_table)

favorable_markers=c('GPC3','PIGY','TMEM150B','MUC13','TREM2','SLC51B','MELK','DTL','TGM3','GPR158')
favorable_markers_short=c('GPC3','TGM3','PIGY','MELK','MUC13')

#reference markers.
reference_genes=read.table('D:/Dropbox/Work/bispecific_markers_project/scrna/PanglaoDB_markers_27_Mar_2020.tsv',sep='\t',header=T)
table(reference_genes$organ)

#Surface_markers

###############################################################
###############################################################
###############################################################
#TUMOR
###############################################################
###############################################################
###############################################################
setwd('D:/Dropbox/Work/bispecific_markers_project/scrna/tumor/')
setwd('~/Dropbox/Work/bispecific_markers_project/scrna/tumor/')


data=read.table('D:/Dropbox/Work/scRNA/tumor_data_all_norm.txt')
data=read.table('~/Dropbox/Work/scRNA/tumor_data_all_norm.txt')
#example_data=data[,sample(1:34414,5000)]
# Load the tumor_seurat dataset

# Initialize the Seurat object with the raw (non-normalized data).
tumor_seurat <- CreateSeuratObject(counts = data, project = "hcc10k", min.cells = 3, min.features = 200)
rm(data)
tumor_seurat[["percent.mt"]] <- PercentageFeatureSet(tumor_seurat, pattern = "^MT-")
VlnPlot(tumor_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


tumor_seurat <- subset(tumor_seurat, subset = nFeature_RNA > 200 & percent.mt < 5)
tumor_seurat <- NormalizeData(tumor_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
tumor_seurat <- FindVariableFeatures(tumor_seurat, selection.method = "vst", nfeatures = 2000)


tumor_classification <- SingleR(test = as.matrix(log2(tumor_seurat@assays[["RNA"]]@counts)),
                                ref = HumanPrimaryCellAtlasData(), labels = HumanPrimaryCellAtlasData()$label.main)

tumor_seurat$labels=tumor_classification$labels
table(tumor_classification$labels)


tumor_seurat <- subset(tumor_seurat, subset = nFeature_RNA > 0 & percent.mt < 50 )
tumor_seurat <- NormalizeData(tumor_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
tumor_seurat <- FindVariableFeatures(tumor_seurat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(tumor_seurat)
tumor_seurat <- ScaleData(tumor_seurat, features = all.genes)
tumor_seurat <- FindVariableFeatures(tumor_seurat)
tumor_seurat <- RunPCA(tumor_seurat, features = VariableFeatures(object = tumor_seurat))
#print(all_data_seurat[["pca"]], dims = 1:5, nfeatures = 5)
#VizDimLoadings(tumor_seurat, dims = 1:2, reduction = "pca")
#heps <- JackStraw(heps, num.replicate = 100)
#heps <- ScoreJackStraw(heps, dims = 1:50)
tumor_seurat <- FindNeighbors(tumor_seurat, dims = 1:30,verbose=TRUE)

set.seed(123)
tumor_seurat <- FindClusters(tumor_seurat, resolution = 0.8)
head(Idents(tumor_seurat), 5)
tumor_seurat <- RunUMAP(tumor_seurat, dims = 1:30)
ElbowPlot(tumor_seurat,reduction='umap')

pdf('liver_tumor.pdf',width=14)
DimPlot(tumor_seurat, reduction = "umap",label=T,group.by = 'labels')+ggtitle('UMAP of HCC dataset')
#heps=subset(heps,idents!=c('8','16','11','15'))

saveRDS(tumor_seurat,file='tumor_seurat.Rds')

#CombinePlots(plots = list(
#  FeaturePlot(healthy_liver, features = 'LINC01983',label=T), 
#  FeaturePlot(healthy_liver, features = 'CRABP1',label=T), 
#  FeaturePlot(healthy_liver, features = 'APOC1',label=T) 
#))
#Glomerulus=4


favorable_markers=c('GPC3','PIGY','TMEM150B','MUC13','TREM2','SLC51B','MELK','DTL','TGM3','GPR158')
favorable_markers_short=c('GPC3','TGM3','PIGY','MELK','MUC13')

#FeaturePlot(tumor_seurat, features = favorable_markers_short)
#FeaturePlot(tumor_seurat, features = favorable_markers)
#FeaturePlot(tumor_seurat, features = 'MGC14156')
#dev.off()


#pdf('liver_tumor_vlnplots.pdf')
VlnPlot(tumor_seurat, features = favorable_markers[favorable_markers!='PIGY'],combine=F,group.by = 'labels')
#VlnPlot(tumor_seurat, features = 'GPC3',combine=F,group.by = 'labels',remove.legend=TRUE)+ theme(legend.position = 'none')
#VlnPlot(kidney_scrna, features = favorable_markers_short,ncol=2)
dev.off()


row.names(data)[grepl('HPM',row.names(data))]
row.names(data)[grepl('MGC',row.names(data))]
row.names(data)[grepl('PIG',row.names(data))]
row.names(data)[grepl('PIG',row.names(data))]
row.names(data)[grepl('PIG',row.names(data))]





###############################################################
###############################################################
###############################################################
#BRAIN
###############################################################
###############################################################
###############################################################
setwd('D:/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/BRAIN/')
setwd('~/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/BRAIN/')


anno=read.csv('GSE138852_covariates.csv',header=T)
anno=subset(anno,oupSample.batchCond=='ct')
row.names(anno)=anno$X
anno[anno$oupSample.cellType=='astro','cell.name']='Astrocytes'
anno[anno$oupSample.cellType=='doublet','cell.name']="Unknown"
anno[anno$oupSample.cellType=='endo','cell.name']='Endothelial'
anno[anno$oupSample.cellType=='mg','cell.name']='Microglia'
anno[anno$oupSample.cellType=='oligo','cell.name']='Oligodendrocytes'
anno[anno$oupSample.cellType=='OPC','cell.name']="Oligodendrocyte progenitor"
anno[anno$oupSample.cellType=='unID','cell.name']="Unknown"



head(anno)
counts=read.csv('GSE138852_counts.csv')
row.names(counts)=counts$X
counts$X=NULL
counts=counts[row.names(anno)]


healthy_brain <- CreateSeuratObject(counts = counts, project = "hcc10k", min.cells = 3, min.features = 10,meta.data = anno)
healthy_brain[["percent.mt"]] <- PercentageFeatureSet(healthy_brain, pattern = "^MT-")
VlnPlot(healthy_brain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 0)
healthy_brain <- subset(healthy_brain, subset = nFeature_RNA > 0 & percent.mt < 50 )
healthy_brain <- NormalizeData(healthy_brain, normalization.method = "LogNormalize", scale.factor = 10000)
healthy_brain <- FindVariableFeatures(healthy_brain, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(healthy_brain)
healthy_brain <- ScaleData(healthy_brain, features = all.genes)
healthy_brain <- FindVariableFeatures(healthy_brain)
healthy_brain <- RunPCA(healthy_brain, features = VariableFeatures(object = healthy_brain))
#print(all_data_seurat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(healthy_brain, dims = 1:2, reduction = "pca")
#heps <- JackStraw(heps, num.replicate = 100)
#heps <- ScoreJackStraw(heps, dims = 1:50)
healthy_brain <- FindNeighbors(healthy_brain, dims = 1:30,verbose=TRUE)

set.seed(123)
healthy_brain <- FindClusters(healthy_brain, resolution = 0.8)
head(Idents(healthy_brain), 5)
healthy_brain <- RunUMAP(healthy_brain, dims = 1:30)
ElbowPlot(healthy_brain,reduction='umap')



pdf('brain_vital.pdf')
DimPlot(healthy_brain, reduction = "umap",label=T)+ggtitle('UMAP of healthy brain tissues')
DimPlot(healthy_brain, reduction = "umap",label=T,group.by = 'cell.name')+ggtitle('UMAP of healthy brain tissues')
#heps=subset(heps,idents!=c('8','16','11','15'))

#Glomerulus=4


favorable_markers=c('GPC3','PIGY','TMEM150B','MUC13','TREM2','SLC51B','MELK','DTL','TGM3','GPR158')
#favorable_markers_short=c('GPC3','TGM3','PIGY','MELK','MUC13')
favorable_markers=favorable_markers[favorable_markers%in%row.names(counts)]


#FeaturePlot(healthy_lung, features = favorable_markers_short)
FeaturePlot(healthy_brain, features = favorable_markers)

VlnPlot(healthy_brain, features = favorable_markers,combine=F,group.by = 'cell.name')
#VlnPlot(kidney_scrna, features = favorable_markers_short,ncol=2)
dev.off()
#################################
#habib17.processed.h5ad
##############################
habib17.processed.h5ad


######################################
#brain habib
#########################################

healthy_brain=readH5AD('habib17.processed.h5ad')
brain_meta=data.frame(cellid=colnames(healthy_brain),celltype=healthy_brain$CellType)
row.names(brain_meta)=brain_meta$cellid
healthy_brain <- CreateSeuratObject(counts=as.matrix(assay(healthy_brain)),meta=brain_meta)
#FeaturePlot(healthy_brain, features = favorable_markers)
healthy_brain[["percent.mt"]] <- PercentageFeatureSet(healthy_brain, pattern = "^MT-")
VlnPlot(healthy_brain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 0)
healthy_brain <- subset(healthy_brain, subset = nFeature_RNA > 0 & percent.mt < 50 )
healthy_brain <- NormalizeData(healthy_brain, normalization.method = "LogNormalize", scale.factor = 10000)
healthy_brain <- FindVariableFeatures(healthy_brain, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(healthy_brain)
healthy_brain <- ScaleData(healthy_brain, features = all.genes)
healthy_brain <- FindVariableFeatures(healthy_brain)
healthy_brain <- RunPCA(healthy_brain, features = VariableFeatures(object = healthy_brain))
#print(all_data_seurat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(healthy_brain, dims = 1:2, reduction = "pca")
#heps <- JackStraw(heps, num.replicate = 100)
#heps <- ScoreJackStraw(heps, dims = 1:50)
healthy_brain <- FindNeighbors(healthy_brain, dims = 1:30,verbose=TRUE)

set.seed(123)
healthy_brain <- FindClusters(healthy_brain, resolution = 0.8)
head(Idents(healthy_brain), 5)
healthy_brain <- RunUMAP(healthy_brain, dims = 1:30)
ElbowPlot(healthy_brain,reduction='umap')

FeaturePlot(healthy_brain, features = favorable_markers)

VlnPlot(healthy_brain, features = favorable_markers,combine=F,group.by = 'cellType')

DimPlot(healthy_brain, reduction = "umap",label=T)+ggtitle('UMAP of healthy brain tissues')
###############################################################
###############################################################
###############################################################
#HEART
###############################################################
###############################################################
###############################################################
setwd('D:/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/HEART/')
setwd('~/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/HEART/')

counts=read.csv('GSE109816_normal_heart_umi_matrix.csv')
row.names(counts)=counts$X
counts$X=NULL
anno=read.table('GSE109816_normal_heart_cell_cluster_info.txt',header=T)
table(anno$CellType)

anno[anno$CellType=='CM','cell.name']='Cardiomyocytes'
anno[anno$CellType=='EC','cell.name']='Endothelial'
anno[anno$CellType=='FB','cell.name']='Fibroblasts'
anno[anno$CellType=='MP','cell.name']='Macrophages'
anno[anno$CellType=='SMC','cell.name']="Smooth muscle"

row.names(anno)=anno$ID

healthy_heart <- CreateSeuratObject(counts = counts, project = "hcc10k", min.cells = 3, min.features = 10,meta.data = anno)
healthy_heart[["percent.mt"]] <- PercentageFeatureSet(healthy_heart, pattern = "^MT-")
VlnPlot(healthy_heart, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 0)
healthy_heart <- subset(healthy_heart, subset = nFeature_RNA > 0 & percent.mt < 50 )
healthy_heart <- NormalizeData(healthy_heart, normalization.method = "LogNormalize", scale.factor = 10000)
healthy_heart <- FindVariableFeatures(healthy_heart, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(healthy_heart)
healthy_heart <- ScaleData(healthy_heart, features = all.genes)
healthy_heart <- FindVariableFeatures(healthy_heart)
healthy_heart <- RunPCA(healthy_heart, features = VariableFeatures(object = healthy_heart))
#print(all_data_seurat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(healthy_heart, dims = 1:2, reduction = "pca")
#heps <- JackStraw(heps, num.replicate = 100)
#heps <- ScoreJackStraw(heps, dims = 1:50)
healthy_heart <- FindNeighbors(healthy_heart, dims = 1:30,verbose=TRUE)

set.seed(123)
healthy_heart <- FindClusters(healthy_heart, resolution = 0.8)
head(Idents(healthy_heart), 5)
healthy_heart <- RunUMAP(healthy_heart, dims = 1:30)
ElbowPlot(healthy_heart,reduction='umap')



pdf('heart_vital.pdf')
DimPlot(healthy_heart, reduction = "umap",label=T)+ggtitle('UMAP of healthy heart tissues')
DimPlot(healthy_heart, reduction = "umap",label=T,group.by = 'cell.name')+ggtitle('UMAP of healthy lung tissues')
#heps=subset(heps,idents!=c('8','16','11','15'))

#Glomerulus=4


favorable_markers=c('GPC3','PIGY','TMEM150B','MUC13','TREM2','SLC51B','MELK','DTL','TGM3','GPR158')
#favorable_markers_short=c('GPC3','TGM3','PIGY','MELK','MUC13')
favorable_markers=favorable_markers[favorable_markers%in%row.names(counts)]


#FeaturePlot(healthy_lung, features = favorable_markers_short)
FeaturePlot(healthy_heart, features = favorable_markers)

VlnPlot(healthy_heart, features = favorable_markers,combine=F,group.by = 'cell.name')
#VlnPlot(kidney_scrna, features = favorable_markers_short,ncol=2)
dev.off()




###############################################################
###############################################################
###############################################################
#LUNG
###############################################################
###############################################################
###############################################################
setwd('D:/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/LUNG/')
setwd('~/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/LUNG/')




load('GSE130148_raw_counts.RData')

cell_clusters=read.table('GSE130148_barcodes_cell_types.txt',header=T,sep='\t')
row.names(cell_clusters)=cell_clusters$cell.barcode

healthy_lung <- CreateSeuratObject(counts = raw_counts, project = "hcc10k", min.cells = 3, min.features = 10,meta.data = cell_clusters)
healthy_lung[["percent.mt"]] <- PercentageFeatureSet(healthy_lung, pattern = "^MT-")
VlnPlot(healthy_lung, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 0)
healthy_lung <- subset(healthy_lung, subset = nFeature_RNA > 0 & percent.mt < 50 )
healthy_lung <- NormalizeData(healthy_lung, normalization.method = "LogNormalize", scale.factor = 10000)
healthy_lung <- FindVariableFeatures(healthy_lung, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(healthy_lung)
healthy_lung <- ScaleData(healthy_lung, features = all.genes)
healthy_lung <- FindVariableFeatures(healthy_lung)
healthy_lung <- RunPCA(healthy_lung, features = VariableFeatures(object = healthy_lung))
#print(all_data_seurat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(healthy_lung, dims = 1:2, reduction = "pca")
#heps <- JackStraw(heps, num.replicate = 100)
#heps <- ScoreJackStraw(heps, dims = 1:50)
healthy_lung <- FindNeighbors(healthy_lung, dims = 1:30,verbose=TRUE)

set.seed(123)
healthy_lung <- FindClusters(healthy_lung, resolution = 0.8)
head(Idents(healthy_lung), 5)
healthy_lung <- RunUMAP(healthy_lung, dims = 1:30)
ElbowPlot(healthy_lung,reduction='umap')



pdf('lung_vital.pdf')
DimPlot(healthy_lung, reduction = "umap",label=T)+ggtitle('UMAP of healthy lung tissues')
DimPlot(healthy_lung, reduction = "umap",label=T,group.by = 'celltype')+ggtitle('UMAP of healthy lung tissues')
#heps=subset(heps,idents!=c('8','16','11','15'))

#Glomerulus=4


favorable_markers=c('GPC3','PIGY','TMEM150B','MUC13','TREM2','SLC51B','MELK','DTL','TGM3','GPR158')
#favorable_markers_short=c('GPC3','TGM3','PIGY','MELK','MUC13')
favorable_markers=favorable_markers[favorable_markers%in%row.names(raw_counts)]


#FeaturePlot(healthy_lung, features = favorable_markers_short)
FeaturePlot(healthy_lung, features = favorable_markers)

VlnPlot(healthy_lung, features = 'HOPX',combine=F,group.by = 'celltype')
#VlnPlot(kidney_scrna, features = favorable_markers_short,ncol=2)
dev.off()



######################################
#Vieira19
#########################################
setwd('D:/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/LUNG/')
setwd('~/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/LUNG/')

lung=readH5AD('vieira19_Alveoli_and_parenchyma_anonymised.processed.h5ad')


###############################################################
###############################################################
###############################################################
#LIVER
###############################################################
###############################################################
###############################################################
setwd('D:/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/LIVER/')
setwd('~/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/LIVER/')

healthy_liver=read.table('D:/Dropbox/Work/scRNA/GSE149614/GSE124395_Normalhumanlivercellatlasdata.txt')
liver_clusters=read.table('D:/Dropbox/Work/scRNA/GSE149614/GSE124395_clusterpartition.txt')
healthy_liver=read.table('~/Dropbox/Work/scRNA/GSE149614/GSE124395_Normalhumanlivercellatlasdata.txt')
liver_clusters=read.table('~/Dropbox/Work/scRNA/GSE149614/GSE124395_clusterpartition.txt')

liver_clusters$cell.type='Unknown'
liver_clusters[liver_clusters$sct.cpart%in%c(5,1,19,3,28,12,18),'cell.type']='NK, NKT, T cells'
liver_clusters[liver_clusters$sct.cpart%in%c(20,9,13,32,23),'cell.type']='Liver sinusoidal end'
liver_clusters[liver_clusters$sct.cpart%in%c(29,33,35,26),'cell.type']='Macrovascular end'
liver_clusters[liver_clusters$sct.cpart%in%c(36,21,9),'cell.type']='Unknown'
liver_clusters[liver_clusters$sct.cpart%in%c(11,17,30,14,8),'cell.type']='Hepatocytes'
liver_clusters[liver_clusters$sct.cpart%in%c(22,38,16,34,37),'cell.type']='B cells'
liver_clusters[liver_clusters$sct.cpart%in%c(25,31,6,23),'cell.type']='Kupfer cells'
liver_clusters[liver_clusters$sct.cpart%in%c(4,7,24,39,11),'cell.type']='EPCAM+'

table(liver_clusters$cell.type)
unique(subset(liver_clusters,cell.type=='A')$sct.cpart)

healthy_liver <- CreateSeuratObject(counts = healthy_liver, project = "hcc10k", min.cells = 0, min.features = 0,meta.data = liver_clusters)
healthy_liver[["percent.mt"]] <- PercentageFeatureSet(healthy_liver, pattern = "^MT-")
VlnPlot(healthy_liver, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 0)
healthy_liver <- subset(healthy_liver, subset = nFeature_RNA > 0 & percent.mt < 50 )
healthy_liver <- NormalizeData(healthy_liver, normalization.method = "LogNormalize", scale.factor = 10000)
healthy_liver <- FindVariableFeatures(healthy_liver, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(healthy_liver)
healthy_liver <- ScaleData(healthy_liver, features = all.genes)
healthy_liver <- FindVariableFeatures(healthy_liver)
healthy_liver <- RunPCA(healthy_liver, features = VariableFeatures(object = healthy_liver))
#print(all_data_seurat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(healthy_liver, dims = 1:2, reduction = "pca")
#heps <- JackStraw(heps, num.replicate = 100)
#heps <- ScoreJackStraw(heps, dims = 1:50)
healthy_liver <- FindNeighbors(healthy_liver, dims = 1:30,verbose=TRUE)

set.seed(123)
healthy_liver <- FindClusters(healthy_liver, resolution = 0.8)
head(Idents(healthy_liver), 5)
healthy_liver <- RunUMAP(healthy_liver, dims = 1:30)
ElbowPlot(healthy_liver,reduction='umap')

saveRDS(healthy_liver,file='healthy_liver.rds')

pdf('liver_vital.pdf')
CombinePlots(plots = list(
DimPlot(healthy_liver, reduction = "umap",label=T,group.by = 'cell.type')+ggtitle('UMAP of healthy liver tissues'),
DimPlot(healthy_liver, reduction = "umap",label=T)+ggtitle('UMAP of healthy liver tissues')))
#heps=subset(heps,idents!=c('8','16','11','15'))




CombinePlots(plots = list(
  FeaturePlot(healthy_liver, features = 'LINC01983',label=T), 
  FeaturePlot(healthy_liver, features = 'CRABP1',label=T), 
  FeaturePlot(healthy_liver, features = 'APOC1',label=T) 
))
#Glomerulus=4


favorable_markers=c('GPC3','PIGY','TMEM150B','MUC13','TREM2','SLC51B','MELK','DTL','TGM3','GPR158')
#favorable_markers_short=c('GPC3','TGM3','PIGY','MELK','MUC13')

FeaturePlot(healthy_liver, features = favorable_markers_short)
FeaturePlot(healthy_liver, features = favorable_markers)

VlnPlot(healthy_liver, features = favorable_markers,combine=F,group.by = 'cell.type')
#VlnPlot(kidney_scrna, features = favorable_markers_short,ncol=2)
dev.off()







###############################################################
###############################################################
###############################################################
#KIDNEY
###############################################################
###############################################################
###############################################################
setwd('D:/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/KIDNEY/')
setwd('~/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/KIDNEY/')

#kidney_reference=subset(reference_genes,organ=='Kidney'&sensitivity_human>0&specificity_human>0&canonical.marker==1)
#table(kidney_reference$cell.type)
kidney_reference=read.csv('kidney_ref.csv')
kidney_reference=subset(kidney_reference,pct.1>0.6&pct.2<0.2)
table(kidney_reference$cluster)
head(kidney_reference)



counts2=read.csv('~/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/KIDNEY/GSM5224978_counts_A1_S1.csv')
#counts3=read.csv('~/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/KIDNEY/KIDNEY_GSE127136_IgAnephropathy_counts.csv')
counts2=read.csv('D:/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/KIDNEY/GSM5224978_counts_A1_S1.csv')
#counts3=read.csv('D:/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/KIDNEY/KIDNEY_GSE127136_IgAnephropathy_counts.csv')


counts2=counts2[-c(1:2),]
row.names(counts2)=counts2$X
counts2$X=NULL

kidney_scrna <- CreateSeuratObject(counts = counts2, project = "hcc10k", min.cells = 0, min.features = 0)
kidney_scrna[["percent.mt"]] <- PercentageFeatureSet(kidney_scrna, pattern = "^MT-")
VlnPlot(kidney_scrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 0)
kidney_scrna <- subset(kidney_scrna, subset = nFeature_RNA > 0 & percent.mt < 50 )
kidney_scrna <- NormalizeData(kidney_scrna, normalization.method = "LogNormalize", scale.factor = 10000)
kidney_scrna <- FindVariableFeatures(kidney_scrna, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(kidney_scrna)
kidney_scrna <- ScaleData(kidney_scrna, features = all.genes)
kidney_scrna <- FindVariableFeatures(kidney_scrna)
kidney_scrna <- RunPCA(kidney_scrna, features = VariableFeatures(object = kidney_scrna))
#print(all_data_seurat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(kidney_scrna, dims = 1:2, reduction = "pca")
#heps <- JackStraw(heps, num.replicate = 100)
#heps <- ScoreJackStraw(heps, dims = 1:50)
kidney_scrna <- FindNeighbors(kidney_scrna, dims = 1:30,verbose=TRUE)



set.seed(123)
kidney_scrna <- FindClusters(kidney_scrna, resolution = 0.8)
head(Idents(kidney_scrna), 5)
kidney_scrna <- RunUMAP(kidney_scrna, dims = 1:30)
ElbowPlot(kidney_scrna,reduction='umap')


pdf('kindney.pdf')
DimPlot(kidney_scrna, reduction = "umap",label=T)+ggtitle('UMAP of healthy kidney tissues')
#heps=subset(heps,idents!=c('8','16','11','15'))


head(kidney_reference)


kidney_reference=read.csv('kidney_ref.csv')
#kidney_reference=subset(kidney_reference,pct.1>0.6&pct.2<0.2)
table(kidney_reference$cluster)
subset(kidney_reference,cluster==names(table(kidney_reference$cluster)[5]))


DimPlot(kidney_scrna, reduction = "umap",label=T)+ggtitle('UMAP of healthy kidney tissues')
##############################define clusters

CombinePlots(plots = list(
#CD=5
FeaturePlot(kidney_scrna, features = 'LINC01983',label=T), 
#DCT=6
FeaturePlot(kidney_scrna, features = 'CRABP1',label=T), 
#DCT-CNT=3
FeaturePlot(kidney_scrna, features = 'APOC1',label=T) 
))
#Glomerulus=4
CombinePlots(plots = list(
FeaturePlot(kidney_scrna, features = 'NPHS1',label=T),
FeaturePlot(kidney_scrna, features = 'TYRO3',label=T),
FeaturePlot(kidney_scrna, features = 'ATP10A',label=T)
))
#Interstitium =8
CombinePlots(plots = list(
FeaturePlot(kidney_scrna, features = 'MYH11',label=T),
FeaturePlot(kidney_scrna, features = 'RERGL',label=T),
FeaturePlot(kidney_scrna, features = 'MCAM',label=T)
))
#PT pure =2
CombinePlots(plots = list(
FeaturePlot(kidney_scrna, features = 'LINC01874',label=T),
FeaturePlot(kidney_scrna, features = 'ALB',label=T),
FeaturePlot(kidney_scrna, features = 'MCAM',label=T)
))
#TAL pure=7,0
CombinePlots(plots = list(
FeaturePlot(kidney_scrna, features = 'ANKRD2',label=T),
FeaturePlot(kidney_scrna, features = 'TDGF1',label=T),
FeaturePlot(kidney_scrna, features = 'PLAU',label=T),
FeaturePlot(kidney_scrna, features = 'GNG7',label=T),
FeaturePlot(kidney_scrna, features = 'PRR15L',label=T),
FeaturePlot(kidney_scrna, features = 'HSPB7',label=T)
))

new.cluster.ids <- c("TAL pure","TAL/PT","PT pure","DCT-CNT","Glomerulus","CD","DCT","TAL pure","Interstitium")
names(new.cluster.ids) <- levels(kidney_scrna)
kidney_scrna <- RenameIdents(kidney_scrna, new.cluster.ids)
DimPlot(kidney_scrna, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()+ggtitle('Labeled healthy cell populations')
kidney_scrna$cell_type=kidney_scrna@active.ident

saveRDS(kidney_scrna,file='kidney_scrna.rds')

favorable_markers=c('GPC3','PIGY','TMEM150B','MUC13','TREM2','SLC51B','MELK','DTL','TGM3','GPR158')
favorable_markers_short=c('GPC3','TGM3','PIGY','MELK','MUC13')

FeaturePlot(kidney_scrna, features = favorable_markers_short,label=T)
FeaturePlot(kidney_scrna, features = favorable_markers,label=T)

VlnPlot(kidney_scrna, features = favorable_markers,combine=F)
#VlnPlot(kidney_scrna, features = favorable_markers_short,ncol=2)
dev.off()
###############################################################
###############################################################
###############################################################



########prepare for report
##################################################


tumor_HEPS=readRDS('~/Dropbox/Work/bispecific_markers_project/scrna/tumor/filtered_and_merged_heps.rds')
meta=tumor_HEPS@meta.data
heps_from_tumor_counts_boolean=data.matrix(as.data.frame(tumor_HEPS@assays[["RNA"]]
                                               [1:nrow(tumor_HEPS@assays[["RNA"]]),
                                                 1:ncol(tumor_HEPS@assays[["RNA"]])]>0)) #plug to transform counts into boolean matrix




pair_of_interest=c('MUC13','GPC3')


colnames(heps_from_tumor_counts_boolean)=meta$type
heps_from_tumor_counts_boolean=heps_from_tumor_counts_boolean[pair_of_interest,]
heps_from_tumor_counts_boolean=heps_from_tumor_counts_boolean[1,]+heps_from_tumor_counts_boolean[2,]
heps_from_tumor_counts_boolean=table(number=heps_from_tumor_counts_boolean,row.names=meta$type)
heps_from_tumor_counts_boolean=as.data.frame.matrix(heps_from_tumor_counts_boolean)
#colnames(heps_from_tumor_counts_boolean)=paste('Tumor_heps_',colnames(heps_from_tumor_counts_boolean))
colnames(heps_from_tumor_counts_boolean)[1:2]=c('Malignant_HEPs','Non-malignant_HEPs')
rm(tumor_HEPS)
#row.names(subset(meta,type=='case'))
#tumor=tumor_counts_boolean[genes_of_interest[1],row.names(subset(meta,type=='case'))]
#tumor plus healthy from tumor
#tumor_only_counts_boolean=heps_from_tumor_counts_boolean[pair_of_interest,row.names(subset(meta,type=='case'))]
#healthy_from_tumor_counts_boolean=heps_from_tumor_counts_boolean[,row.names(subset(meta,type=='control'))]
#tumor_co_ex

#co_expression_tumor=tumor_only_counts_boolean[1,]+tumor_only_counts_boolean[2,row.names(subset(meta,type=='case'))]
#co_expression_tumor=table(co_expression_tumor)/length(row.names(subset(meta,type=='case')))
#healthy_from_tumor  

#co_expression_healthy_heps_from_tumor=healthy_from_tumor_counts_boolean[1,]+healthy_from_tumor_counts_boolean[2,row.names(subset(meta,type=='control'))]
#co_expression_healthy_heps_from_tumor=table(co_expression_healthy_heps_from_tumor)/length(row.names(subset(meta,type=='control')))
#if(length(co_expression_healthy_heps_from_tumor)<3){
#  co_expression_healthy_heps_from_tumor['2']=0
#}
###########################################
#kidney
kidney_scrna=readRDS('~/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/KIDNEY/kidney_scrna.rds')
kidney_meta=kidney_scrna@meta.data


kidney_boolean=data.matrix(as.data.frame(kidney_scrna@assays[["RNA"]]
                                                         [1:nrow(kidney_scrna@assays[["RNA"]]),
                                                           1:ncol(kidney_scrna@assays[["RNA"]])]>0)) #plug to transform counts into boolean matrix
colnames(kidney_boolean)=kidney_meta$cell_type
kidney_boolean=kidney_boolean[pair_of_interest,]
kidney_boolean=kidney_boolean[1,]+kidney_boolean[2,]
kidney_boolean=table(number=kidney_boolean,row.names=kidney_meta$cell_type)
kidney_boolean=as.data.frame.matrix(kidney_boolean)
colnames(kidney_boolean)=paste('Kidney_',colnames(kidney_boolean))
rm(kidney_scrna)

###########################################
#liver
liver_boolean=readRDS('~/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/LIVER/healthy_liver.rds')
liver_meta=liver_boolean@meta.data


liver_boolean=data.matrix(as.data.frame(liver_boolean@assays[["RNA"]]
                                         [1:nrow(liver_boolean@assays[["RNA"]]),
                                           1:ncol(liver_boolean@assays[["RNA"]])]>0)) #plug to transform counts into boolean matrix
colnames(liver_boolean)=liver_meta$cell.type
liver_boolean=liver_boolean[pair_of_interest,]
liver_boolean=liver_boolean[1,]+liver_boolean[2,]
liver_boolean=table(number=liver_boolean,row.names=liver_meta$cell.type)
liver_boolean=as.data.frame.matrix(liver_boolean)
if(nrow(liver_boolean)<3){
  liver_boolean['2',]=0
}
colnames(liver_boolean)=paste('Healthy_liver_',colnames(liver_boolean),sep='')
rm(liver_scrna)
###########################################
#lungs


######overall table

for_heatmap=cbind(heps_from_tumor_counts_boolean,kidney_boolean,liver_boolean)

#col_sum=apply(for_heatmap,2,sum)
for_heatmap=apply(for_heatmap,2,function(x){x/sum(x)})
row.names(for_heatmap)=c('None','Single positive',paste0(pair_of_interest[1],'+',pair_of_interest[2],'+'))
par(mar=c(7,4,4,2))
heatmap.2(as.matrix(for_heatmap),col='bluered',dendrogram = 'none',cexCol=0.5,cexRow=0.7,las=2,srtCol=45,srtRow=45)

#################################################