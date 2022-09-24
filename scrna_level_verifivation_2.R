library(dplyr)
library(Seurat)
library(ggrepel)
library(ggplot2)

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

########################
setwd('~/Dropbox/Work/bispecific_markers_project/scrna/')
setwd('D:/Dropbox/Work/bispecific_markers_project/scrna/')
setwd('tumor')
cell_list=list.files(pattern='.Rds')
setwd('..')

cell_list=c('BM.Rds',
'B_cell.Rds',
'Chondrocytes.Rds',
'CMP.Rds',
'DC.Rds',
'Epithelial_cells.Rds',
'Fibroblasts.Rds',
'GMP.Rds',
'Hepatocytes.Rds',
'HSC_CD34+.Rds',
'Macrophage.Rds',
'MEP.Rds',
'Monocyte.Rds',
'Neutrophils.Rds',
'NK_cell.Rds',
'Osteoblasts.Rds',
'Pre-B_cell_CD34-.Rds',
'Pro-B_cell_CD34+.Rds',
'Smooth_muscle_cells.Rds',
'Tissue_stem_cells.Rds',
'T_cells.Rds')

i=cell_list[15] #macro
i=cell_list[19] #NK
i=cell_list[24] #T


i
a=0
for(i in cell_list){
load(paste0('tumor/',i))
tumor=temp
tumor$pheno='tumor'
load(paste0('healthy/',i))
healthy=temp
healthy$pheno='healthy'
temp=merge(tumor,healthy)
temp=merge(tumor,healthy)
if(sum(table(healthy$orig.ident))>50&sum(table(tumor$orig.ident))>50){
temp_DE=FindMarkers(temp, ident.1='tumor',group.by='pheno',logfc.threshold=0) 
temp_DE$gene=row.names(temp_DE)
temp_DE_surface=temp_DE[row.names(temp_DE)%in%unique(c(result_table$antigen_1,result_table$antigen_2)),]
temp_DE_surface=na.omit(temp_DE_surface)
write.table(temp_DE_surface,file=paste0('results/',i,'_DE_surface.txt'))
write.table(temp_DE,file=paste0('results/',i,'_DE.txt'))

all.genes <- rownames(temp)
temp <- ScaleData(temp, features = all.genes)
temp <- FindVariableFeatures(temp)
temp <- RunPCA(temp, features = VariableFeatures(object = temp),npcs=min(50,ncol(temp)-1))
#print(all_data_seurat[["pca"]], dims = 1:5, nfeatures = 5)
#VizDimLoadings(temp, dims = 1:2, reduction = "pca")
#heps <- JackStraw(heps, num.replicate = 100)
#heps <- ScoreJackStraw(heps, dims = 1:50)
temp <- FindNeighbors(temp, dims = 1:length(temp$pca))
temp <- FindClusters(temp, resolution = 1)
temp <- RunUMAP(temp, dims = 1:(length(temp$pca)-1),n.neighbors=(length(temp$pca)-1))

temp_DE$diffexpressed='NO'
temp_DE[temp_DE$p_val_adj<0.05&temp_DE$avg_log2FC>0.5,'diffexpressed']='UP'
temp_DE[temp_DE$p_val_adj<0.05&temp_DE$avg_log2FC<(-0.5),'diffexpressed']='DOWN'
temp_DE[row.names(temp_DE_surface),'diffexpressed']='Surface'
#label DE genes
temp_DE$delabel <- NA
#temp_DE$delabel[temp_DE$diffexpressed != "NO"] <- temp_DE$gene[temp_DE$diffexpressed != "NO"]
temp_DE$delabel[temp_DE$diffexpressed == "Surface"&temp_DE$diffexpressed != "NO"] <- temp_DE$gene[temp_DE$diffexpressed == "Surface"&temp_DE$diffexpressed != "NO"]
#fancy volcano plot

pdf(paste0('results/',i,'_plots.pdf'))
for(x in row.names(temp_DE_surface)){
print(FeaturePlot(temp, features = x,label=T)+ ggtitle(paste0('Expression of ',x,' in ',gsub('.Rds','',i))))
}
print(DimPlot(temp, reduction = "umap",label=F,group.by='pheno')+ ggtitle(gsub('.Rds','',i)))
print(DimPlot(temp, reduction = "umap",label=F,group.by='orig.ident')+ ggtitle(gsub('.Rds','',i)))
print(ggplot(data=temp_DE, aes(x=avg_log2FC , y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black","darkgreen", "red")) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red"))

dev.off()
}
a=a+1
cat(round(a/length(cell_list),2),'\n')
}






#fancy volcano plot
ggplot(data=DE_toy, aes(x=avg_log2FC , y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")




















immune_cells=c("B_cell.Rds","Macrophage.Rds","NK_cell.Rds","Monocyte.Rds","T_cells.Rds")

load(file=cell_list[12])
current_type=temp
current_type$pheno='control'
current_type$pheno[grepl('T',current_type$orig.ident)]='case'
table(current_type$pheno)


current_type_DE=FindMarkers(current_type, ident.1='case',group.by='pheno',logfc.threshold=0) #
current_type_DE$Gene_id=row.names(DE_toy)


#label DE genes
DE_toy$delabel <- NA
DE_toy$delabel[DE_toy$diffexpressed != "NO"] <- DE_toy$hgnc_symbol[DE_toy$diffexpressed != "NO"]

#fancy volcano plot
ggplot(data=DE_toy, aes(x=avg_log2FC , y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")


