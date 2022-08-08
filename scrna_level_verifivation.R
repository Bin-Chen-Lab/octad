
#primary cells
setwd('D:/Dropbox/Work/bispecific_markers_project/scrna/primary')
setwd('~/Dropbox/Work/bispecific_markers_project/scrna/primary')

result_table=read.table('D:/Dropbox/Work/bispecific_markers_project/bulk/test_SVM_new_4.txt')
result_table=read.table('~/Dropbox/Work/bispecific_markers_project/bulk/test_SVM_new_4.txt')

result_table=na.omit(result_table)
plot(density(result_table$pair_score),lwd=3,col='red',main='Distribution of pair score')
abline(v=quantile(result_table$pair_score,.90),lwd=4,col='blue')
result_table=subset(result_table,pair_score>quantile(result_table$pair_score,.90))
result_table=subset(result_table,case_greater=='TRUE_TRUE')
result_table=result_table[order(result_table$pair_score,decreasing = T),]
head(result_table)

plot(density(result_table$pair_score),lwd=3,col='red',main='Distribution of pair score')

result_table=na.omit(result_table)

malignant_heps=readRDS('filtered_and_merged_heps.rds')
malignant_heps$Malign.type[is.na(malignant_heps$Malign.type)]='healthy'
table(malignant_heps$Malign.type)
malignant_heps=magic(malignant_heps,genes='all_genes')
malignant_heps@active.assay <- 'MAGIC_RNA'



CombinePlots(plots = list(
  #  FeaturePlot(heps, features = c('CD4'),label=T),
  #  FeaturePlot(heps, features = c('AMBP'),label=T),
  FeaturePlot(malignant_heps, features = c('GPC3'),label=T),
#  FeaturePlot(malignant_heps, features = c('AMBP'),label=T),
#  DimPlot(malignant_heps, reduction = "umap",label=T,group.by='type'),
#  DimPlot(heps, reduction = "umap",label=T),
  DimPlot(malignant_heps, reduction = "umap",label=T,group.by='Malign.type')
)
)






load('healthy_HEPS.Rds')
# Initialize the Seurat object with the raw (non-normalized data).
healthy_heps <- CreateSeuratObject(counts = healthy_heps_counts, project = "hcc10k", min.cells = 3, min.features = 200)
healthy_heps[["percent.mt"]] <- PercentageFeatureSet(healthy_heps, pattern = "^MT-")
VlnPlot(healthy_heps, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


healthy_heps <- subset(healthy_heps, subset = nFeature_RNA > 200 & percent.mt < 5)
healthy_heps <- NormalizeData(healthy_heps, normalization.method = "LogNormalize", scale.factor = 10000)
healthy_heps <- FindVariableFeatures(healthy_heps, selection.method = "vst", nfeatures = 2000)



merged_heps=merge(malignant_heps,healthy_heps)
rm(healthy_heps)
rm(healthy_heps_counts)
rm(mmalignant_heps)


















heps_DE=FindMarkers(heps, ident.1=names(heps$type[heps$type=='case']), ident.2=names(heps$type[heps$type=='control']))
heps_DE=heps_DE[order(heps_DE$avg_log2FC,decreasing=T),]
#heps_DE_filtered=heps_DE[unique(c(result_table$antigen_1,result_table$antigen_2))%in%row.names(heps_DE),]
heps_DE_filtered=heps_DE[unique(c(result_table$antigen_1,result_table$antigen_2)),]
heps_DE_filtered=na.omit(heps_DE_filtered)

FeaturePlot(object = heps, features = c('MUC13', 'GPC3'), cols = c("grey", "red", "blue"), blend=TRUE,ncol=2,combine=T)
FeaturePlot(object = heps, features = c('MUC13', 'GPC3'), cols = c("grey", "red", "blue"), blend=TRUE,ncol=2,combine=T)
FeaturePlot(object = heps, features = c('MUC13', 'GPC3'), cols = c("grey", "red", "blue"), blend=TRUE,ncol=2,combine=T)

CombinePlots(plots = list(
  #  FeaturePlot(heps, features = c('CD4'),label=T),
  #  FeaturePlot(heps, features = c('AMBP'),label=T),
  FeaturePlot(heps, features = c('GPC3'),label=T),
  FeaturePlot(heps, features = c('AMBP'),label=T),
  DimPlot(heps, reduction = "umap",label=T,group.by='type'),
  DimPlot(heps, reduction = "umap",label=T),
  DimPlot(heps, reduction = "umap",label=T,group.by='orig.ident')
)
)



CombinePlots(plots = list(
  #  FeaturePlot(heps, features = c('CD4'),label=T),
  #  FeaturePlot(heps, features = c('AMBP'),label=T),
  DimPlot(heps, reduction = "umap",label=T),
  DimPlot(heps, reduction = "umap",label=T,group.by='orig.ident')
)
)
FindMarkers(heps, ident.1=names(heps$type[heps$type=='case']), ident.2=names(heps$type[heps$type=='control']))
#result_table=subset(result_table,angle_cos>0.7)
plot(density(result_table$pair_score))


result_table$antigen_1=gsub('-','_',result_table$antigen_1)
result_table$antigen_2=gsub('-','_',result_table$antigen_2)
#colnames(train)=gsub('-','_',colnames(train))
#immune cells
setwd('~/Dropbox/Work/bispecific_markers_project/scrna/immune')
load('tumor_T_cells_counts.Rds')
load('healthy_T_cells_counts.Rds')
#bridge method from pseudobulk

t_cells=cbind(healthy_T_cells_counts,tumor_T_cells_counts)
meta=data.frame(sample=c(colnames(healthy_T_cells_counts),colnames(tumor_T_cells_counts)),
                pheno=c(rep('Healthy_T',ncol(healthy_T_cells_counts)),rep('Tumor_T',ncol(tumor_T_cells_counts))))
row.names(meta)=meta$sample


  t_cells <- CreateSeuratObject(counts = t_cells, project = "hcc10k", min.cells = 3, min.features = 200)
  t_cells[["percent.mt"]] <- PercentageFeatureSet(t_cells, pattern = "^MT-")
VlnPlot(t_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



t_cells <- subset(t_cells, subset = nFeature_RNA > 200 & percent.mt < 5)
t_cells <- NormalizeData(t_cells, normalization.method = "LogNormalize", scale.factor = 10000)
t_cells <- FindVariableFeatures(t_cells, selection.method = "vst", nfeatures = 2000)


tumor_classification <- SingleR(test = as.matrix(log2(t_cells@assays[["RNA"]]@counts)),
                                ref = HumanPrimaryCellAtlasData(), labels = HumanPrimaryCellAtlasData()$label.main)


table(t_cells$labels)
t_cells$pruned.labels=tumor_classification$pruned.labels

t_cells=subset(t_cells, subset = pruned.labels == "T_cells"&orig.ident!='HCC04T')
t_cells=subset(t_cells, subset = pruned.labels == "T_cells"&orig.ident!='HCC04N')


all.genes <- rownames(t_cells )
t_cells  <- ScaleData(t_cells)
t_cells  <- FindVariableFeatures(t_cells )
t_cells  <- RunPCA(t_cells , features = VariableFeatures(object = t_cells ))
#print(all_data_seurat[["pca"]], dims = 1:5, nfeatures = 5)
#VizDimLoadings(tumor_seurat , dims = 1:2, reduction = "pca")
#heps <- JackStraw(heps, num.replicate = 100)
#heps <- ScoreJackStraw(heps, dims = 1:50)
#ElbowPlot(heps,reduction='umap')
tumor_seurat  <- FindNeighbors(tumor_seurat , dims = 1:30)
tumor_seurat  <- FindClusters(tumor_seurat , resolution = 1)
head(Idents(tumor_seurat ), 5)
tumor_seurat  <- RunUMAP(tumor_seurat , dims = 1:30)

#heps=subset(heps,idents!=c('8','16','11','15'))

pdf('cd4+gpc3+ambp.pdf',height=12,width=12)
CombinePlots(plots = list(
  #  FeaturePlot(heps, features = c('CD4'),label=T),
  #  FeaturePlot(heps, features = c('AMBP'),label=T),
  FeaturePlot(t_cells, features = c('GPC3'),label=T),
  DimPlot(tumor_seurat, reduction = "umap",label=T,group.by='pruned.labels')
))
dev.off()