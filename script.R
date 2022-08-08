load("/Users/chekali1/Dropbox/Work/scRNA/HCC_final_pipeline_scRNA/tumor_HEPS.Rds")
load("/Users/chekali1/Dropbox/Work/scRNA/HCC_final_pipeline_scRNA/healthy_HEPS.Rds") #unix
setwd('/Users/chekali1/Dropbox/Work/scRNA/malignancy_prediction/')

load("D:/Dropbox/Work/scRNA/HCC_final_pipeline_scRNA/healthy_HEPS.Rds")  #win
load("D:/Dropbox/Work/scRNA/HCC_final_pipeline_scRNA/tumor_HEPS.Rds")
setwd('D:/Dropbox/Work/scRNA/malignancy_prediction/')


merged=cbind(healthy_heps_counts[1:1000],tumor_heps_counts[1:1000])



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


dir.create('test_6')

barcodes=colnames(merged)
write.table(barcodes,file='test_6/filtered_feature_bc_matrix/barcodes.tsv',quote=F,row.names=F,col.names = F,sep='\t')
features=row.names(merged)
write.table(features,file='test_6/filtered_feature_bc_matrix/features.tsv',quote=F,col.names = F,sep='\t')
merged=as.matrix(merged)
#counts=counts[,apply(counts,2,mean)>0.06345945]
merged <- as(as.matrix(merged), "sparseMatrix") 
writeMM(merged,file='test_6/filtered_feature_bc_matrix/matrix.mtx')


healthy_liver=read.table('D:/Dropbox/Work/scRNA/GSE149614/GSE124395_Normalhumanlivercellatlasdata.txt')
liver_clusters=read.table('D:/Dropbox/Work/scRNA/GSE149614/GSE124395_clusterpartition.txt')

healthy_liver=read.table('/Users/chekali1/Dropbox/Work/scRNA/GSE149614/GSE124395_Normalhumanlivercellatlasdata.txt')
liver_clusters=read.table('/Users/chekali1/Dropbox/Work/scRNA/GSE149614/GSE124395_clusterpartition.txt')

liver_clusters$cell_name=row.names(liver_clusters)
liver_clusters=liver_clusters[liver_clusters$sct.cpart%in%c(11,17,14),]
healthy_liver=healthy_liver[liver_clusters$cell_name]


huh7=read.table("D:/Dropbox/Work/scRNA/GSE149614/HuH7.txt")
huh7=read.table("D:/Dropbox/Work/scRNA/GSE149614/HuH7.txt")

huh7=read.table("/Users/chekali1/Dropbox/Work/scRNA/GSE149614/HuH7.txt")
huh7=read.table("/Users/chekali1/Dropbox/Work/scRNA/GSE149614/HuH7.txt")
huh7=huh7[!grepl('3',colnames(huh7))]
#HCC_from_huh7=huh7[grepl('3',colnames(huh7))]
#huh1=huh7[grepl('1',colnames(huh7))]
huh7_for_reference=huh7[grepl('2',colnames(huh7))]
#reference=cbind(huh1[sample(1:ncol(huh1),200)],huh7[sample(1:ncol(huh7),200)],healthy_from_huh7[sample(1:ncol(healthy_from_huh7),200)])
#annotation=data.frame(sample=colnames(reference),phenotype=c(rep('huh1-like',200),rep('huh7-like',200),rep('healthy_heps',200)))


#a=Read10X('test_4/filtered_feature_bc_matrix')
huh7_seurat=CreateSeuratObject(huh7)
#huh7_seurat$Sample=c(rep('huh1',894),rep('huh7',2082),rep('HCC',673))
huh7_seurat$Sample=c(rep('huh1',894),rep('huh7',2082))
huh7_seurat$Batch=1
healthy_liver_seurat=CreateSeuratObject(counts=healthy_heps_counts)
healthy_liver_seurat$Sample='healthy_HEPS'
healthy_liver_seurat$Batch=2
huh7_huh1_hcc_healthy_liver_seurat=merge(huh7_seurat,healthy_liver_seurat)

#Define true classes
#huh7_huh1_hcc_healthy_liver_seurat$True_classified=c(rep('malignant',3649),rep('nonMalignant',2672))
huh7_huh1_hcc_healthy_liver_seurat$True_classified=c(rep('malignant',2976),rep('nonMalignant',2672))
#####################################

huh7_huh1_hcc_healthy_liver_seurat[["percent.mt"]] <- PercentageFeatureSet(huh7_huh1_hcc_healthy_liver_seurat, pattern = "^MT-")
huh7_huh1_hcc_healthy_liver_seurat <- subset(huh7_huh1_hcc_healthy_liver_seurat, subset = nFeature_RNA > 200 & percent.mt < 5)
huh7_huh1_hcc_healthy_liver_seurat <- NormalizeData(huh7_huh1_hcc_healthy_liver_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
huh7_huh1_hcc_healthy_liver_seurat <- FindVariableFeatures(huh7_huh1_hcc_healthy_liver_seurat, selection.method = "vst", nfeatures = 2000)


all.genes <- rownames(huh7_huh1_hcc_healthy_liver_seurat)
huh7_huh1_hcc_healthy_liver_seurat <- ScaleData(huh7_huh1_hcc_healthy_liver_seurat, features = all.genes)
huh7_huh1_hcc_healthy_liver_seurat <- FindVariableFeatures(huh7_huh1_hcc_healthy_liver_seurat)
huh7_huh1_hcc_healthy_liver_seurat <- RunPCA(huh7_huh1_hcc_healthy_liver_seurat, features = VariableFeatures(object = huh7_huh1_hcc_healthy_liver_seurat))
#print(all_data_seurat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(huh7_huh1_hcc_healthy_liver_seurat, dims = 1:2, reduction = "pca")
#heps <- JackStraw(heps, num.replicate = 100)
#heps <- ScoreJackStraw(heps, dims = 1:50)
ElbowPlot(huh7_huh1_hcc_healthy_liver_seurat,reduction='umap')
huh7_huh1_hcc_healthy_liver_seurat <- FindNeighbors(huh7_huh1_hcc_healthy_liver_seurat, dims = 1:30)
huh7_huh1_hcc_healthy_liver_seurat <- FindClusters(huh7_huh1_hcc_healthy_liver_seurat, resolution = 0.5)
head(Idents(huh7_huh1_hcc_healthy_liver_seurat), 5)
huh7_huh1_hcc_healthy_liver_seurat <- RunUMAP(huh7_huh1_hcc_healthy_liver_seurat, dims = 1:30)
huh7_huh1_hcc_healthy_liver_seurat <- RunPCA(huh7_huh1_hcc_healthy_liver_seurat)
huh7_huh1_hcc_healthy_liver_seurat <- RunHarmony(huh7_huh1_hcc_healthy_liver_seurat, "Batch")
huh7_huh1_hcc_healthy_liver_seurat <- RunUMAP(huh7_huh1_hcc_healthy_liver_seurat, dims = 1:30)




dir.create('test_9')
dir.create('test_9/filtered_feature_bc_matrix')

counts=huh7_huh1_hcc_healthy_liver_seurat@assays[["RNA"]]@counts
barcodes=colnames(counts)
write.table(barcodes,file='test_9/filtered_feature_bc_matrix/barcodes.tsv',quote=F,row.names=F,col.names = F,sep='\t')
features=row.names(counts)
reference_genes=read.table('test_5/filtered_feature_bc_matrix/features.tsv')
reference_genes=reference_genes[reference_genes$V2 %in% features,]
reference_genes=reference_genes[!duplicated(reference_genes$V2),]
write.table(reference_genes,file='test_9/filtered_feature_bc_matrix/features.tsv',quote=F,col.names = F, row.names=F, sep='\t')
counts=counts[reference_genes$V2 %in% rownames(counts),]
counts=counts[rownames(counts) %in% reference_genes$V2,]
counts=as.matrix(counts)
#counts=counts[,apply(counts,2,mean)>0.06345945]
counts <- as(as.matrix(counts), "sparseMatrix") 
writeMM(counts,file='test_9/filtered_feature_bc_matrix/matrix.mtx')


dataPath <- "test_9"     # The path of cell ranger processed data
savePath <- "test_9"  # A path to save the results
sampleName <- "test_9"          # The sample name
authorName <- ""           # The author name to mark the report

# Run scStatistics
stat.results <- runScStatistics(
  dataPath = dataPath,
  savePath = savePath,
  sampleName = sampleName,
  authorName = authorName
)


dataPath <- "test_9"      # The path of cell ranger processed data
statPath <- "test_9"   # The path of the scStatistics results
savePath <- "test_9"   # A path to save the results
sampleName <- "test_9"           # The sample name
authorName <- ""            

# Run scAnnotation
anno.results <- runScAnnotation(
  dataPath = dataPath,
  statPath = statPath,
  savePath = savePath,
  authorName = authorName,
  sampleName = sampleName,
  geneSet.method = "average",cnv.ref.data=healthy_liver, # or "GSVA",
  bool.runExprProgram=F,
  bool.runDoublet=F,
#  bool.runCellClassify=F,
  bool.runCellCycle=F,
  bool.runStemness=F,
  bool.runGeneSets=F,
  bool.runDiffExpr=F,
  bool.runInteraction=F
)

malignancy=read.table('/Users/chekali1/Dropbox/Work/scRNA/malignancy_prediction/test_8/cellAnnotation.txt',header=T)
malignancy=read.table('D:/Dropbox/Work/scRNA/malignancy_prediction/test_9/cellAnnotation.txt',header=T)
malignancy=subset(malignancy,select=c('barcodes','Malign.score','Malign.type'))
row.names(malignancy)=malignancy$barcodes
malignancy=malignancy[colnames(huh7_huh1_hcc_healthy_liver_seurat),]

unique(colnames(huh7_huh1_hcc_healthy_liver_seurat)==malignancy$barcodes)

huh7_huh1_hcc_healthy_liver_seurat$Malign.score=malignancy$Malign.score
huh7_huh1_hcc_healthy_liver_seurat$Malign.type=malignancy$Malign.type
#$huh7_huh1_hcc_healthy_liver_seurat


CombinePlots(plots = list(
DimPlot(huh7_huh1_hcc_healthy_liver_seurat, reduction = "umap",label=T,group.by='Malign.type'),
DimPlot(huh7_huh1_hcc_healthy_liver_seurat, reduction = "umap",label=T,group.by='Sample')
))
#DimPlot(huh7_huh1_hcc_healthy_liver_seurat, reduction = "umap",label=T,group.by='pheno')

##########################################test model performance
df_for_test=data.frame(predicted=as.factor(huh7_huh1_hcc_healthy_liver_seurat$Malign.type),
                                reference=as.factor(huh7_huh1_hcc_healthy_liver_seurat$True_classified))
df_for_test=na.omit(df_for_test_INFERCNV)
sensitivity(df_for_test_INFERCNV$predicted,df_for_test_INFERCNV$reference)
##########################################
sensitivity(as.factor(huh7_huh1_hcc_healthy_liver_seurat$Malign.type), 
            as.factor(huh7_huh1_hcc_healthy_liver_seurat$True_classified))
specificity(as.factor(huh7_huh1_hcc_healthy_liver_seurat$Malign.type), 
            as.factor(huh7_huh1_hcc_healthy_liver_seurat$True_classified))


########################################################################################
#My part
########################################################################################

#Define reference
huh7_for_reference=huh7[grepl('2',colnames(huh7))]
huh7_for_reference=CreateSeuratObject(counts=huh7_for_reference,class='malignant')
huh7_for_reference$Batch=1
  
healthy_liver=read.table('D:/Dropbox/Work/scRNA/GSE149614/GSE124395_Normalhumanlivercellatlasdata.txt')
liver_clusters=read.table('D:/Dropbox/Work/scRNA/GSE149614/GSE124395_clusterpartition.txt')

healthy_liver=read.table('/Users/chekali1/Dropbox/Work/scRNA/GSE149614/GSE124395_Normalhumanlivercellatlasdata.txt')
liver_clusters=read.table('/Users/chekali1/Dropbox/Work/scRNA/GSE149614/GSE124395_clusterpartition.txt')

liver_clusters$cell_name=row.names(liver_clusters)
liver_clusters=liver_clusters[liver_clusters$sct.cpart%in%c(11,17,14),]
healthy_liver=healthy_liver[liver_clusters$cell_name]
healthy_liver=CreateSeuratObject(counts=healthy_liver)
healthy_liver$class='nonMalignant'
healthy_liver$Batch=2

reference_singleR=merge(huh7_for_reference,healthy_liver)

reference_singleR[["percent.mt"]] <- PercentageFeatureSet(reference_singleR, pattern = "^MT-")
reference_singleR <- subset(reference_singleR, subset = nFeature_RNA > 200 & percent.mt < 5)
reference_singleR <- NormalizeData(reference_singleR, normalization.method = "LogNormalize", scale.factor = 10000)
reference_singleR <- FindVariableFeatures(reference_singleR, selection.method = "vst", nfeatures = 2000)


all.genes <- rownames(reference_singleR)
reference_singleR <- ScaleData(reference_singleR, features = all.genes)
reference_singleR <- FindVariableFeatures(reference_singleR)
reference_singleR <- RunPCA(reference_singleR, features = VariableFeatures(object = reference_singleR))
#print(all_data_seurat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(reference_singleR, dims = 1:2, reduction = "pca")
#heps <- JackStraw(heps, num.replicate = 100)
#heps <- ScoreJackStraw(heps, dims = 1:50)
ElbowPlot(reference_singleR,reduction='umap')
reference_singleR <- FindNeighbors(reference_singleR, dims = 1:30)
reference_singleR <- FindClusters(reference_singleR, resolution = 0.5)
head(Idents(reference_singleR), 5)
reference_singleR <- RunUMAP(reference_singleR, dims = 1:30)
reference_singleR <- RunPCA(reference_singleR)
reference_singleR <- RunHarmony(reference_singleR, "Batch")
reference_singleR <- RunUMAP(reference_singleR, dims = 1:30)

reference_singleR$class=c(rep('malignant',2082),rep('nonMalignant',3040))
##########Train classifier:
reference_singleR <- getFeatureSpace(reference_singleR, "class")
reference_singleR <- trainModel(reference_singleR)

#

huh7_huh1_hcc_healthy_liver_seurat <- scPredict(huh7_huh1_hcc_healthy_liver_seurat, reference_singleR)
huh7_huh1_hcc_healthy_liver_seurat$scpred_prediction[huh7_huh1_hcc_healthy_liver_seurat$scpred_prediction=='unassigned']=NA

#######################
#
#######################

CombinePlots(plots = list(
  DimPlot(huh7_huh1_hcc_healthy_liver_seurat, reduction = "umap",label=T,group.by='Malign.type'),
  DimPlot(huh7_huh1_hcc_healthy_liver_seurat, reduction = "umap",label=T,group.by='Sample')
))
#DimPlot(huh7_huh1_hcc_healthy_liver_seurat, reduction = "umap",label=T,group.by='pheno')

##########################################test model performance
df_for_test=data.frame(predicted=as.factor(huh7_huh1_hcc_healthy_liver_seurat$Malign.type),
                       reference=as.factor(huh7_huh1_hcc_healthy_liver_seurat$True_classified))
df_for_test=na.omit(df_for_test_INFERCNV)
sensitivity(df_for_test_INFERCNV$predicted,df_for_test_INFERCNV$reference)
##########################################
sensitivity(as.factor(huh7_huh1_hcc_healthy_liver_seurat$Malign.type), 
            as.factor(huh7_huh1_hcc_healthy_liver_seurat$True_classified))
specificity(as.factor(huh7_huh1_hcc_healthy_liver_seurat$Malign.type), 
            as.factor(huh7_huh1_hcc_healthy_liver_seurat$True_classified))

sensitivity(as.factor(huh7_huh1_hcc_healthy_liver_seurat$scpred_prediction), 
            as.factor(huh7_huh1_hcc_healthy_liver_seurat$True_classified))
specificity(as.factor(huh7_huh1_hcc_healthy_liver_seurat$scpred_prediction), 
            as.factor(huh7_huh1_hcc_healthy_liver_seurat$True_classified))

DimPlot(huh7_huh1_hcc_healthy_liver_seurat, reduction = "umap",label=T,group.by='Sample')
DimPlot(huh7_huh1_hcc_healthy_liver_seurat, reduction = "umap",label=T,group.by='Malign.type')+ggtitle('inferCNV based on liver atlas')
DimPlot(huh7_huh1_hcc_healthy_liver_seurat, reduction = "umap",label=T,group.by='scpred_prediction')+ggtitle('our method based on liver atlas+huh7 cells')
