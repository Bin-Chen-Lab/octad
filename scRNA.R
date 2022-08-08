
library(dplyr)
library(Seurat)
library(garnett)
library(monocle)
library(splitstackshape)
library(celldex)
library(scater)
library(scRNAseq)
library(Matrix)

#classification packages
library(infercnv)
library(SingleR)
library(scPred)

setwd('D:/Dropbox/Work/scRNA/HCC_final_pipeline_scRNA')
setwd('~/Dropbox/Work/scRNA/HCC_final_pipeline_scRNA')
load('tumor_HEPS.Rds')

ну 
load('tumor_T_cells_counts.Rds')
#tumor data

#############################################
#############################################load and classify tumor data
#############################################
data=read.table('../tumor_data_all_norm.txt')

example_data=data[,sample(1:34414,5000)]
# Load the tumor_seurat dataset

# Initialize the Seurat object with the raw (non-normalized data).
tumor_seurat <- CreateSeuratObject(counts = tumor_heps_counts, project = "hcc10k", min.cells = 3, min.features = 200)
tumor_seurat[["percent.mt"]] <- PercentageFeatureSet(tumor_seurat, pattern = "^MT-")
VlnPlot(tumor_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


tumor_seurat <- subset(tumor_seurat, subset = nFeature_RNA > 200 & percent.mt < 5)
tumor_seurat <- NormalizeData(tumor_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
tumor_seurat <- FindVariableFeatures(tumor_seurat, selection.method = "vst", nfeatures = 2000)


tumor_classification <- SingleR(test = as.matrix(log2(tumor_seurat@assays[["RNA"]]@counts)),
                                ref = HumanPrimaryCellAtlasData(), labels = HumanPrimaryCellAtlasData()$label.main)


table(tumor_classification$labels)

#save T cells and
tumor_heps=subset(tumor_classification,labels=='Hepatocytes')
tumor_T_cells=subset(tumor_classification,labels=='T_cells')

tumor_heps_counts=data[row.names(tumor_heps)]
tumor_T_cells_counts=data[row.names(tumor_T_cells)]

save(tumor_heps_counts,file='tumor_HEPS.Rds')
save(tumor_T_cells_counts,file='tumor_T_cells_counts.Rds')

load('tumor_HEPS.Rds')
load('tumor_T_cells_counts.Rds')

#############################################
#############################################load and classify healthy data
#############################################
data=read.table('../healthy_data_norm.txt')
# Load the tumor_seurat dataset

# Initialize the Seurat object with the raw (non-normalized data).
healthy_seurat <- CreateSeuratObject(counts = data, project = "hcc10k", min.cells = 3, min.features = 200)
healthy_seurat[["percent.mt"]] <- PercentageFeatureSet(healthy_seurat, pattern = "^MT-")
VlnPlot(healthy_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


healthy_seurat <- subset(healthy_seurat, subset = nFeature_RNA > 200 & percent.mt < 5)
healthy_seurat <- NormalizeData(healthy_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
healthy_seurat <- FindVariableFeatures(healthy_seurat, selection.method = "vst", nfeatures = 2000)


healthy_classification <- SingleR(test = as.matrix(log2(healthy_seurat@assays[["RNA"]]@counts)), 
                                  ref = HumanPrimaryCellAtlasData(), labels = HumanPrimaryCellAtlasData()$label.main)


table(healthy_classification$labels)



table(t_cells$labels)
healthy_seurat$pruned.labels=healthy_classification$pruned.labels

t_cells=subset(t_cells, subset = pruned.labels == "T_cells"&orig.ident!='HCC04T')
t_cells=subset(t_cells, subset = pruned.labels == "T_cells"&orig.ident!='HCC04N')




#save T cells and
for(i in names(table(healthy_classification$labels))){
  temp=subset(healthy_seurat, subset = pruned.labels == i&orig.ident!='HCC04N')
  save(temp,file=paste0(i,'.Rds'))
  
}
healthy_heps=subset(healthy_classification,labels=='Hepatocytes')
healthy_T_cells=subset(healthy_classification,labels=='T_cells')

healthy_heps_counts=data[row.names(healthy_heps)]
healthy_T_cells_counts=data[row.names(healthy_T_cells)]

save(healthy_heps_counts,file='healthy_HEPS.Rds')
save(healthy_T_cells_counts,file='healthy_T_cells_counts.Rds')

load('healthy_HEPS.Rds')
load('healthy_T_cells_counts.Rds')
#############################################
#############################################load and classify malignant cells
#############################################

load('healthy_HEPS.Rds')
load('tumor_HEPS.Rds')

data=cbind(healthy_heps_counts,tumor_heps_counts)
meta=data.frame(sample=c(colnames(healthy_heps_counts),colnames(tumor_heps_counts)),
                pheno=c(rep('Healthy_derived_HEPs',ncol(healthy_heps_counts)),rep('Tumor_derived_HEPs',ncol(tumor_heps_counts))))
row.names(meta)=meta$sample

heps <- CreateSeuratObject(counts = data, project = "hcc10k", min.cells = 3, min.features = 200,meta.data = meta)
heps[["percent.mt"]] <- PercentageFeatureSet(heps, pattern = "^MT-")
VlnPlot(heps, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


heps <- subset(heps, subset = nFeature_RNA > 200 & percent.mt < 5)
heps <- NormalizeData(heps, normalization.method = "LogNormalize", scale.factor = 10000)
heps <- FindVariableFeatures(heps, selection.method = "vst", nfeatures = 2000)


all.genes <- rownames(heps)
heps <- ScaleData(heps, features = all.genes)
heps <- FindVariableFeatures(heps)
heps <- RunPCA(heps, features = VariableFeatures(object = heps))
#print(all_data_seurat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(heps, dims = 1:2, reduction = "pca")
#heps <- JackStraw(heps, num.replicate = 100)
#heps <- ScoreJackStraw(heps, dims = 1:50)
ElbowPlot(heps,reduction='umap')
heps <- FindNeighbors(heps, dims = 1:30)
heps <- FindClusters(heps, resolution = 1)
head(Idents(heps), 5)
heps <- RunUMAP(heps, dims = 1:30)

#heps=subset(heps,idents!=c('8','16','11','15'))

pdf('cd4+gpc3+ambp.pdf',height=12,width=12)
CombinePlots(plots = list(
#  FeaturePlot(heps, features = c('CD4'),label=T),
#  FeaturePlot(heps, features = c('AMBP'),label=T),
  FeaturePlot(heps, features = c('GPC3'),label=T),
  DimPlot(heps, reduction = "umap",label=T,group.by='pheno')
))
dev.off()

pdf('gpc3_ambp_coexpression.pdf',height=5,width=20)
FeaturePlot(object = heps, features = c('AMBP', 'GPC3'), cols = c("grey", "red", "blue"), blend=TRUE,ncol=2,combine=T)

dev.off()

