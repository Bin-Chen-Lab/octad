library(Matrix)
library(knitr)
library(Seurat)
library(ggplot2)
library(scCancer)
library(harmony)
library(dplyr) 
library(garnett)
library(monocle)
library(splitstackshape)
library(celldex)
library(scater)
library(scRNAseq)
library(Matrix)
library(caret)
#classification packages
library(infercnv)
library(SingleR)
library(scPred)


load("/Users/chekali1/Dropbox/Work/scRNA/HCC_final_pipeline_scRNA/tumor_HEPS.Rds")
load("/Users/chekali1/Dropbox/Work/scRNA/HCC_final_pipeline_scRNA/healthy_HEPS.Rds") #unix
setwd('/Users/chekali1/Dropbox/Work/scRNA/malignancy_prediction/')

load("D:/Dropbox/Work/scRNA/HCC_final_pipeline_scRNA/healthy_HEPS.Rds")  #win
load("D:/Dropbox/Work/scRNA/HCC_final_pipeline_scRNA/tumor_HEPS.Rds")
setwd('D:/Dropbox/Work/scRNA/malignancy_prediction/')


healthy_liver=read.table('D:/Dropbox/Work/scRNA/GSE149614/GSE124395_Normalhumanlivercellatlasdata.txt')
liver_clusters=read.table('D:/Dropbox/Work/scRNA/GSE149614/GSE124395_clusterpartition.txt')

healthy_liver=read.table('/Users/chekali1/Dropbox/Work/scRNA/GSE149614/GSE124395_Normalhumanlivercellatlasdata.txt')
liver_clusters=read.table('/Users/chekali1/Dropbox/Work/scRNA/GSE149614/GSE124395_clusterpartition.txt')

liver_clusters$cell_name=row.names(liver_clusters)
liver_clusters=liver_clusters[liver_clusters$sct.cpart%in%c(11,17,14),]
healthy_liver=healthy_liver[liver_clusters$cell_name]



counts=cbind(healthy_heps_counts,tumor_heps_counts)

dir.create('test_11')
dir.create('test_11/filtered_feature_bc_matrix')

#counts=cbind(healthy_heps_counts,tumor_heps_counts)
counts=tumor_heps_counts
barcodes=colnames(counts)
write.table(barcodes,file='test_11/filtered_feature_bc_matrix/barcodes.tsv',quote=F,row.names=F,col.names = F,sep='\t')
features=row.names(counts)
reference_genes=read.table('test_5/filtered_feature_bc_matrix/features.tsv')
reference_genes=reference_genes[reference_genes$V2 %in% features,]
reference_genes=reference_genes[!duplicated(reference_genes$V2),]
write.table(reference_genes,file='test_11/filtered_feature_bc_matrix/features.tsv',quote=F,col.names = F, row.names=F, sep='\t')
counts=counts[reference_genes$V2 %in% rownames(counts),]
counts=counts[rownames(counts) %in% reference_genes$V2,]
counts=as.matrix(counts)
#counts=counts[,apply(counts,2,mean)>0.06345945]
counts <- as(as.matrix(counts), "sparseMatrix") 
writeMM(counts,file='test_11/filtered_feature_bc_matrix/matrix.mtx')


dataPath <- "test_11"     # The path of cell ranger processed data
savePath <- "test_11"  # A path to save the results
sampleName <- "test_11"          # The sample name
authorName <- ""           # The author name to mark the report

# Run scStatistics
stat.results <- runScStatistics(
  dataPath = dataPath,
  savePath = savePath,
  sampleName = sampleName,
  authorName = authorName
)


dataPath <- "test_11"      # The path of cell ranger processed data
statPath <- "test_11"   # The path of the scStatistics results
savePath <- "test_11"   # A path to save the results
sampleName <- "test_11"           # The sample name
authorName <- ""            

# Run scAnnotation
anno.results <- runScAnnotation(
  dataPath = dataPath,
  statPath = statPath,
  savePath = savePath,
  authorName = authorName,
  sampleName = sampleName,
  geneSet.method = "average",cnv.ref.data=healthy_heps_counts, # or "GSVA",
  bool.runExprProgram=F,
  bool.runDoublet=F,
  #  bool.runCellClassify=F,
  bool.runCellCycle=F,
  bool.runStemness=F,
  bool.runGeneSets=F,
  bool.runDiffExpr=F,
  bool.runInteraction=F
)

malignancy=read.table('/Users/chekali1/Dropbox/Work/scRNA/malignancy_prediction/test_11/cellAnnotation.txt',header=T)
malignancy=read.table('D:/Dropbox/Work/scRNA/malignancy_prediction/test_11/cellAnnotation.txt',header=T)
malignancy=subset(malignancy,select=c('barcodes','Malign.score','Malign.type'))
row.names(malignancy)=malignancy$barcodes
malignancy=malignancy[colnames(tumor_heps_counts),]


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

tumor_seurat$class=tumor_classification$pruned.labels
tumor_seurat$Malign.type=malignancy$Malign.type
tumor_seurat$Malign.score=malignancy$Malign.score


all.genes <- rownames(tumor_seurat)
tumor_seurat <- ScaleData(tumor_seurat, features = all.genes)
tumor_seurat <- FindVariableFeatures(tumor_seurat)
tumor_seurat <- RunPCA(tumor_seurat, features = VariableFeatures(object = tumor_seurat))
#print(all_data_seurat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(tumor_seurat, dims = 1:2, reduction = "pca")
#heps <- JackStraw(heps, num.replicate = 100)
#heps <- ScoreJackStraw(heps, dims = 1:50)
ElbowPlot(tumor_seurat,reduction='umap')
tumor_seurat <- FindNeighbors(tumor_seurat, dims = 1:30)
tumor_seurat <- FindClusters(tumor_seurat, resolution = 1)
head(Idents(tumor_seurat), 5)
tumor_seurat <- RunUMAP(tumor_seurat, dims = 1:30)

tumor_seurat_malignant_heps=subset(tumor_seurat, subset = orig.ident!='HCC04T'&orig.ident!='HCC04N')



#tumor_seurat_malignant_heps=subset(tumor_seurat, subset = Malign.type=='malignant'&orig.ident!='HCC04T')

CombinePlots(plots = list(
  #  FeaturePlot(heps, features = c('CD4'),label=T),
  #  FeaturePlot(heps, features = c('AMBP'),label=T),
  FeaturePlot(tumor_seurat_malignant_heps, features = c('GPC3'),label=T),
  FeaturePlot(tumor_seurat_malignant_heps, features = c('AMBP'),label=T),
  DimPlot(tumor_seurat_malignant_heps, reduction = "umap",label=T,group.by='Malign.type'),
  DimPlot(tumor_seurat_malignant_heps, reduction = "umap",label=T,group.by='orig.ident')
))


all.genes <- rownames(tumor_seurat_malignant_heps)
tumor_seurat_malignant_heps <- ScaleData(tumor_seurat_malignant_heps, features = all.genes)
tumor_seurat_malignant_heps <- FindVariableFeatures(tumor_seurat_malignant_heps)
tumor_seurat_malignant_heps <- RunPCA(tumor_seurat_malignant_heps, features = VariableFeatures(object = tumor_seurat_malignant_heps))
#print(all_data_seurat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(tumor_seurat_malignant_heps, dims = 1:2, reduction = "pca")
#heps <- JackStraw(heps, num.replicate = 100)
#heps <- ScoreJackStraw(heps, dims = 1:50)
ElbowPlot(tumor_seurat_malignant_heps,reduction='umap')

tumor_seurat_malignant_heps <- RunHarmony(tumor_seurat_malignant_heps, "orig.ident")
tumor_seurat_malignant_heps <- RunUMAP(tumor_seurat_malignant_heps, reduction = "harmony")

tumor_seurat_malignant_heps <- FindNeighbors(tumor_seurat_malignant_heps, dims = 1:30)
tumor_seurat_malignant_heps <- FindClusters(tumor_seurat_malignant_heps, resolution = 1)
head(Idents(tumor_seurat_malignant_heps), 5)
tumor_seurat_malignant_heps <- RunUMAP(tumor_seurat_malignant_heps, dims = 1:30)

tumor_seurat_malignant_heps$type='case'

###############################################################

# Initialize the Seurat object with the raw (non-normalized data).
healthy_seurat <- CreateSeuratObject(counts = healthy_heps_counts, project = "hcc10k", min.cells = 3, min.features = 200)
healthy_seurat[["percent.mt"]] <- PercentageFeatureSet(healthy_seurat, pattern = "^MT-")
VlnPlot(healthy_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


healthy_seurat <- subset(healthy_seurat, subset = nFeature_RNA > 200 & percent.mt < 5)
healthy_seurat <- NormalizeData(healthy_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
healthy_seurat <- FindVariableFeatures(healthy_seurat, selection.method = "vst", nfeatures = 2000)


healthy_classification <- SingleR(test = as.matrix(log2(healthy_seurat@assays[["RNA"]]@counts)),
                                ref = HumanPrimaryCellAtlasData(), labels = HumanPrimaryCellAtlasData()$label.main)


table(healthy_classification$labels)

healthy_seurat$class=healthy_classification$pruned.labels


all.genes <- rownames(healthy_seurat)
healthy_seurat <- ScaleData(healthy_seurat, features = all.genes)
healthy_seurat <- FindVariableFeatures(healthy_seurat)
healthy_seurat <- RunPCA(healthy_seurat, features = VariableFeatures(object = healthy_seurat))
#print(all_data_seurat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(healthy_seurat, dims = 1:2, reduction = "pca")
#heps <- JackStraw(heps, num.replicate = 100)
#heps <- ScoreJackStraw(heps, dims = 1:50)
ElbowPlot(healthy_seurat,reduction='umap')
healthy_seurat <- FindNeighbors(healthy_seurat, dims = 1:30)
healthy_seurat <- FindClusters(healthy_seurat, resolution = 1)
head(Idents(healthy_seurat), 5)
healthy_seurat <- RunUMAP(healthy_seurat, dims = 1:30)


healthy_seurat=subset(healthy_seurat, subset = class == "Hepatocytes"&orig.ident!='HCC04N')

CombinePlots(plots = list(
  #  FeaturePlot(heps, features = c('CD4'),label=T),
  #  FeaturePlot(heps, features = c('AMBP'),label=T),
  FeaturePlot(healthy_seurat, features = c('GPC3'),label=T),
  FeaturePlot(healthy_seurat, features = c('AMBP'),label=T),
  DimPlot(healthy_seurat, reduction = "umap",label=T,group.by='orig.ident')
))
healthy_seurat$type='control'


merged_heps=merge(tumor_seurat_malignant_heps,healthy_seurat)
all.genes <- rownames(merged_heps)
merged_heps <- ScaleData(merged_heps, features = all.genes)
merged_heps <- FindVariableFeatures(merged_heps)
merged_heps <- RunPCA(merged_heps, features = VariableFeatures(object = merged_heps))
#print(all_data_seurat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(merged_heps, dims = 1:2, reduction = "pca")
#heps <- JackStraw(heps, num.replicate = 100)
#heps <- ScoreJackStraw(heps, dims = 1:50)
ElbowPlot(merged_heps,reduction='umap')
merged_heps <- FindNeighbors(merged_heps, dims = 1:30)
merged_heps <- FindClusters(merged_heps, resolution = 1)
head(Idents(merged_heps), 5)
merged_heps <- RunUMAP(merged_heps, dims = 1:30)

merged_heps <- RunHarmony(merged_heps, "orig.ident",max.iter.harmony = 20,
                          max.iter.cluster = 30,assay.use="RNA")
merged_heps <- RunUMAP(merged_heps, reduction = "harmony",dims=1:30)

CombinePlots(plots = list(
  #  FeaturePlot(heps, features = c('CD4'),label=T),
  #  FeaturePlot(heps, features = c('AMBP'),label=T),
  FeaturePlot(merged_heps, features = c('GPC3'),label=T),
  FeaturePlot(merged_heps, features = c('AMBP'),label=T),
  DimPlot(merged_heps, reduction = "harmony",label=T,group.by='orig.ident'),
  DimPlot(merged_heps, reduction = "harmony",label=T,group.by='type')
))


saveRDS(merged_heps,file='filtered_and_merged_heps.rds')
