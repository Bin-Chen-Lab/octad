# Load libraries
library(scater)
library(Seurat)
library(SingleCellExperiment)


library(tidyverse) 
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)

library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)

# Extract raw counts and metadata to create SingleCellExperiment object
counts <- temp@assays$RNA@counts 

metadata <- temp@meta.data

# Set up metadata as desired for aggregation and DE analysis
metadata$seurat_clusters  <- factor(temp@active.ident)
metadata$orig.ident  <- factor(metadata$orig.ident)
metadata$Malign.type   <- factor(metadata$Malign.type )
metadata$pheno   <- factor(metadata$pheno )

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

summed <- aggregateAcrossCells(sce, 
                               id=colData(sce)[,c("seurat_clusters", "orig.ident",'pheno')])
summed

# Creating up a DGEList object for use in edgeR:
library(edgeR)
y <- DGEList(counts(summed), samples=colData(summed))
y

discarded <- summed$ncells < 10
y <- y[,!discarded]
summary(discarded)
keep <- filterByExpr(y, group=y$samples$Malign.type.1)
y <- y[keep,]
summary(keep)
y <- calcNormFactors(y)
y$samples

colnames(y$counts)=paste(y$samples$orig.ident.1,y$samples$pheno.1,y$samples$seurat_clusters.1,sep='_')
y$counts[1:10,1:10]

save(y,file=paste0(i,'pseudobulk.rds'))



####################################################
load('T_cells.Rdspseudobulk.rds')
t_cells=y



t_cells$samples$pheno.1
t_cells$samples$color='#F8766D'
t_cells$samples$color[t_cells$samples$pheno.1=='tumor']='#00C19A'
pca=prcomp(t_cells$counts,center=T,scale=T)
pca=as.data.frame(pca$x)
plot(pca$PC1,pca$PC2,pch=21,bg=t_cells$samples$color,main='T cells pseudobulk',xlab='PC1',ylab='PC2')

load('Macrophage.Rdspseudobulk.rds')
macro=y
macro$samples$pheno.1
macro$samples$color='#F8766D'
macro$samples$color[macro$samples$pheno.1=='tumor']='#00C19A'
pca=prcomp(macro$counts)
pca=as.data.frame(pca$x)
plot(pca$PC1,pca$PC2,pch=21,bg=t_cells$samples$color,main=paste0(i),xlab='PC1',ylab='PC2')


load('heps_pseudobulk.RDS')
heps=y

heps$samples$Malign.type.1
heps$samples$color='#F8766D'
heps$samples$color[heps$samples$Malign.type.1=='malignant']='#00C19A'
pca=prcomp(heps$counts,center=T,scale=T)
pca=as.data.frame(pca$x)
plot(pca$PC1,pca$PC2,pch=21,bg=heps$samples$color,main='HEP cells pseudobulk',xlab='PC1',ylab='PC2')


heps$counts=log2(heps$counts)
plot(heps$counts['PLVAP',],heps$counts['GPC3',],pch=21,bg=heps$samples$color)















load('~/Dropbox/Work/bispecific_markers_project/scrna/immune/healthy_T_cells_counts.Rds')
load('~/Dropbox/Work/bispecific_markers_project/scrna/immune/tumor_T_cells_counts.Rds')

T_cells=cbind(healthy_T_cells_counts,tumor_T_cells_counts)
metadata=data.frame(sample=colnames(T_cells),orig.ident=gsub("\\_.*",'',colnames(T_cells)))
row.names(metadata)=metadata$sample

T_cells <- CreateSeuratObject(counts = T_cells, project = "hcc10k", min.cells = 3, min.features = 200,meta.data = metadata)
T_cells[["percent.mt"]] <- PercentageFeatureSet(T_cells, pattern = "^MT-")
VlnPlot(T_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


T_cells <- subset(T_cells, subset = nFeature_RNA > 200 & percent.mt < 5)
T_cells <- NormalizeData(T_cells, normalization.method = "LogNormalize", scale.factor = 10000)
T_cells <- FindVariableFeatures(T_cells, selection.method = "vst", nfeatures = 2000)


all.genes <- rownames(T_cells)
T_cells <- ScaleData(T_cells, features = all.genes)
T_cells <- FindVariableFeatures(T_cells)
T_cells <- RunPCA(T_cells, features = VariableFeatures(object = T_cells))
#print(all_data_seurat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(T_cells, dims = 1:2, reduction = "pca")
#heps <- JackStraw(heps, num.replicate = 100)
#heps <- ScoreJackStraw(heps, dims = 1:50)
ElbowPlot(T_cells,reduction='umap')
T_cells <- FindNeighbors(T_cells, dims = 1:30)
T_cells <- FindClusters(T_cells, resolution = 1)
head(Idents(T_cells), 5)
T_cells <- RunUMAP(T_cells, dims = 1:30)


# Extract raw counts and metadata to create SingleCellExperiment object
counts <- T_cells@assays$RNA@counts 


metadata <- T_cells@meta.data

# Set up metadata as desired for aggregation and DE analysis
metadata$seurat_clusters  <- factor(T_cells@active.ident)
metadata$orig.ident  <- factor(T_cells$orig.ident)
metadata$Malign.type='healthy'
metadata$Malign.type[grepl('T',metadata$orig.ident.1)]='malignant'
# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

summed <- aggregateAcrossCells(sce, 
                               id=colData(sce)[,c("seurat_clusters", "orig.ident",'Malign.type')],average=T)
summed

# Creating up a DGEList object for use in edgeR:
library(edgeR)
x <- DGEList(counts(summed), samples=colData(summed))
x
x$samples$Malign.type.1[grepl('T',x$samples$orig.ident.1)]='tumor'
discarded <- summed$ncells < 10
x <- x[,!discarded]
summary(discarded)
keep <- filterByExpr(x, group=x$samples$Malign.type.1)
x <- x[keep,]
summary(keep)
x<- calcNormFactors(x)
x$samples

colnames(x$counts)=paste(x$samples$orig.ident.1,x$samples$Malign.type.1,x$samples$seurat_clusters.1,sep='_')
x$counts[1:10,1:10]

test_T=x$counts[ grep('CD3',row.names(x$counts)),1:10]
test_T=data.frame(CD3G=x$counts['CD3D',],state=as.factor(x$samples$Malign.type.1))
colnames(test_hep)=y$samples$Malign.type.1


test_hep=data.frame(GPC3=y$counts['GPC3',],state=y$samples$Malign.type.1)



boxplot(test_T$CD3G~test_T$state)
boxplot(test_hep$GPC3~test_hep$state)
#colnames(test_hep)=y$samples$Malign.type.1

# Extract raw counts and metadata to create SingleCellExperiment object
counts <- tumor_seurat_malignant_heps@assays$RNA@counts 


metadata <- tumor_seurat_malignant_heps@meta.data

# Set up metadata as desired for aggregation and DE analysis
metadata$seurat_clusters  <- factor(heps@active.ident)
metadata$orig.ident  <- factor(metadata$orig.ident)
metadata$pheno   <- factor(metadata$pheno )

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

# Identify groups for aggregation of counts
groups <- colData(sce)[, c("seurat_clusters", "orig.ident",'pheno')]

# Explore the raw counts for the dataset

## Check the assays present
assays(sce)

## Explore the raw counts for the dataset
dim(counts(sce))

counts(sce)[1:6, 1:6]

## Explore the cellular metadata for the dataset
dim(colData(sce))

head(colData(sce))

# Named vector of cluster names
kids <- purrr::set_names(levels(sce$seurat_clusters))
kids

# Total number of clusters
nk <- length(kids)
nk

# Named vector of sample names
sids <- purrr::set_names(levels(sce$orig.ident))

# Total number of samples 
ns <- length(sids)
ns


# Generate sample level metadata

## Determine the number of cells per sample
table(sce$orig.ident)

## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$orig.ident))

## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$orig.ident)

## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  select(-"seurat_clusters")
ei


# Perform QC if not already performed
dim(sce)

# Calculate quality control (QC) metrics
sce <- perCellQCMetrics(sce)

# Get cells w/ few/many detected genes
sce$is_outlier <- isOutlier(
  metric = sce$total,
  nmads = 2, type = "both", log = TRUE)

# Remove outlier cells
sce=subset(sce,is_outlier==FALSE)
#sce <- sce[, !sce$is_outlier]
dim(sce)

## Remove lowly expressed genes which have less than 10 cells with any counts
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]

dim(sce)


# Aggregate the counts per sample_id and cluster_id

# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("seurat_clusters", "orig.ident",'pheno')]

# Aggregate across cluster-sample groups
pseudomatrix <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

pseudomatrix=as.data.frame(pseudomatrix)
pseudomatrix$pheno='control'
pseudomatrix[grepl('Tumor_derived_HEPs',row.names(pseudomatrix)),'pheno']='case'
table(pseudomatrix$pheno)




