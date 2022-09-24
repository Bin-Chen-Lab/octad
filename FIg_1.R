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



#load results table
result_table=read.table('D:/Dropbox/Work/bispecific_markers_project/bulk/test_SVM_new_4.txt',header=T)
result_table=read.table('~/Dropbox/Work/bispecific_markers_project/bulk/test_SVM_new_4.txt',header=T)



result_table=na.omit(result_table)
plot(density(result_table$pair_score),lwd=3,col='red',main='Distribution of pair score')
abline(v=quantile(result_table$pair_score,.90),lwd=4,col='blue')
result_table=subset(result_table,pair_score>quantile(result_table$pair_score,.90))
result_table=subset(result_table,case_greater=='TRUE_TRUE')
result_table=result_table[order(result_table$pair_score,decreasing = T),]
head(result_table)

result_table=result_table[1:20,]


setwd('~/Dropbox/Work/bispecific_markers_project/scrna/')
setwd('D:/Dropbox/Work/bispecific_markers_project/scrna/')


cell_classifier=read.csv('~/Dropbox/Work/bispecific_markers_project/scrna/Cell_classification.csv')
cell_classifier=read.csv('D:/Dropbox/Work/bispecific_markers_project/scrna/Cell_classification.csv')
cell_classifier=rbind(cell_classifier,c('Neuronal','Neuronal','Unknown','Brain'))

#prepare data
tumor_HEPS=readRDS('~/Dropbox/Work/bispecific_markers_project/scrna/tumor/filtered_and_merged_heps.rds')
tumor_HEPS=readRDS('D:/Dropbox/Work/bispecific_markers_project/scrna/tumor/filtered_and_merged_heps.rds')
meta=tumor_HEPS@meta.data
heps_from_tumor_counts_boolean_whole=data.matrix(as.data.frame(tumor_HEPS@assays[["RNA"]]
                                                               [1:nrow(tumor_HEPS@assays[["RNA"]]),
                                                                 1:ncol(tumor_HEPS@assays[["RNA"]])]>0)) 
colnames(heps_from_tumor_counts_boolean_whole)=meta$type
rm(tumor_HEPS)

#kidney
kidney_scrna=readRDS('~/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/KIDNEY/kidney_scrna_2.rds')
kidney_scrna=readRDS('D:/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/KIDNEY/kidney_scrna_2.rds')
kidney_meta=kidney_scrna@meta.data
kidney_meta$cell_name=row.names(kidney_meta)
head(kidney_meta)
colnames(kidney_meta)[colnames(kidney_meta)=='annotation']='cell.type'
#idney_meta$Organ='Kidney'
dim(kidney_meta)
kidney_meta=merge(kidney_meta,subset(cell_classifier,Organ=='Kidney'),by.x='cell.type',by.y='Level.3',all=T)
dim(kidney_meta)
row.names(kidney_meta)=kidney_meta$cell_name
kidney_scrna@meta.data=kidney_meta
saveRDS(kidney_scrna,file='~/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/KIDNEY/kidney_scrna_2.rds')



kidney_boolean_whole=data.matrix(as.data.frame(kidney_scrna@assays[["RNA"]]
                                         [1:nrow(kidney_scrna@assays[["RNA"]]),
                                           1:ncol(kidney_scrna@assays[["RNA"]])]>0)) #plug to transform counts into boolean matrix
colnames(kidney_boolean_whole)=kidney_meta$annotation
rm(kidney_scrna)
#liver
liver_boolean_whole=readRDS('~/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/LIVER/healthy_liver.rds')
liver_boolean_whole=readRDS('D:/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/LIVER/healthy_liver.rds')
liver_meta=liver_boolean_whole@meta.data
liver_meta$cell_name=row.names(liver_meta)
head(liver_meta)
#colnames(liver_meta)[colnames(liver_meta)=='annotation']='cell.type'
#liver_meta$Organ='Liver'
dim(liver_meta)
liver_meta=merge(liver_meta,subset(cell_classifier,Organ=='Liver'),by.x='cell.type',by.y='Level.3',all=T)
liver_meta=unique(liver_meta)
dim(liver_meta)
row.names(liver_meta)=liver_meta$cell_name
liver_boolean_whole@meta.data=liver_meta
saveRDS(liver_boolean_whole,file='~/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/LIVER/healthy_liver.rds')


liver_boolean_whole=data.matrix(as.data.frame(liver_boolean_whole@assays[["RNA"]]
                                        [1:nrow(liver_boolean_whole@assays[["RNA"]]),
                                          1:ncol(liver_boolean_whole@assays[["RNA"]])]>0)) #plug to transform counts into boolean matrix
colnames(liver_boolean_whole)=liver_meta$cell.type
###############brain
#brain
brain_boolean_whole=readRDS('~/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/BRAIN/healthy_brain.rds')
brain_boolean_whole=readRDS('D:/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/BRAIN/healthy_brain.rds')
brain_meta=brain_boolean_whole@meta.data
brain_meta$cell_name=row.names(brain_meta)
#brain_meta$Organ='Brain'
dim(brain_meta)
brain_meta=merge(brain_meta,subset(cell_classifier,Organ=='Brain'),by.x='cell.name',by.y='Level.3',all=T)
colnames(brain_meta)[colnames(brain_meta)=='cell.name']='cell.type'
dim(brain_meta)
row.names(brain_meta)=brain_meta$cell_name
brain_boolean_whole@meta.data=brain_meta
saveRDS(brain_boolean_whole,file='~/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/BRAIN/healthy_brain.rds')


brain_boolean_whole=data.matrix(as.data.frame(brain_boolean_whole@assays[["RNA"]]
                                        [1:nrow(brain_boolean_whole@assays[["RNA"]]),
                                          1:ncol(brain_boolean_whole@assays[["RNA"]])]>0)) #plug to transform counts into boolean matrix
colnames(brain_boolean_whole)=brain_meta$cell.name
####################
#lung
lung_boolean_whole=readRDS('~/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/LUNG/healthy_lung.rds')
lung_boolean_whole=readRDS('D:/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/LUNG/healthy_lung.rds')
lung_meta=lung_boolean_whole@meta.data
lung_meta$cell_name=row.names(lung_meta)
head(lung_meta)
#lung_meta$Organ='Lung'
dim(lung_meta)
lung_meta=merge(lung_meta,subset(cell_classifier,Organ=='Lung'),by.x='celltype',by.y='Level.3',all=T)
dim(lung_meta)
colnames(lung_meta)[colnames(lung_meta)=='celltype']='cell.type'
row.names(lung_meta)=lung_meta$cell_name
lung_boolean_whole@meta.data=lung_meta
saveRDS(lung_boolean_whole,file='~/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/LUNG/healthy_lung.rds')


lung_boolean_whole=data.matrix(as.data.frame(lung_boolean_whole@assays[["RNA"]]
                                       [1:nrow(lung_boolean_whole@assays[["RNA"]]),
                                         1:ncol(lung_boolean_whole@assays[["RNA"]])]>0)) #plug to transform counts into boolean matrix
colnames(lung_boolean_whole)=lung_meta$celltype

#################heart
heart_boolean_whole=readRDS('~/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/HEART/healthy_heart.rds')
heart_boolean_whole=readRDS('D:/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/HEART/healthy_heart.rds')
heart_meta=heart_boolean_whole@meta.data
colnames(heart_meta)[colnames(heart_meta)=='cell.name']='cell.type'
heart_meta$cell_name=row.names(heart_meta)
#heart_meta$Organ='Heart'
dim(heart_meta)
heart_meta=merge(heart_meta,subset(cell_classifier,Organ=='Heart'),by.x='cell.type',by.y='Level.3',all=T)
dim(heart_meta)
row.names(heart_meta)=heart_meta$cell_name
heart_boolean_whole@meta.data=heart_meta
saveRDS(heart_boolean_whole,file='~/Dropbox/Work/bispecific_markers_project/scrna/REFERENCE_HEALTHY/HEART/healthy_heart.rds')


heart_boolean_whole=data.matrix(as.data.frame(heart_boolean_whole@assays[["RNA"]]
                                        [1:nrow(heart_boolean_whole@assays[["RNA"]]),
                                          1:ncol(heart_boolean_whole@assays[["RNA"]])]>0)) #plug to transform counts into boolean matrix
colnames(heart_boolean_whole)=heart_meta$cell.name


i=1
for(i in 1:nrow(result_table)){
  
pair_of_interest=c(result_table[i,1],result_table[i,2])

###########################################################################################
################TUMOR
#plug to transform counts into boolean matrix

heps_from_tumor_counts_boolean=heps_from_tumor_counts_boolean_whole[row.names(heps_from_tumor_counts_boolean_whole)%in%pair_of_interest,]
if(is.null(nrow(heps_from_tumor_counts_boolean))){
  heps_from_tumor_counts_boolean=rbind(heps_from_tumor_counts_boolean,0)
}else if(nrow(heps_from_tumor_counts_boolean)<1){
  heps_from_tumor_counts_boolean=heps_from_tumor_counts_boolean_whole[1:2,]
  heps_from_tumor_counts_boolean[1,]=0
  heps_from_tumor_counts_boolean[2,]=0
}


heps_from_tumor_counts_boolean=heps_from_tumor_counts_boolean[1,]+heps_from_tumor_counts_boolean[2,]
heps_from_tumor_counts_boolean=table(number=heps_from_tumor_counts_boolean,row.names=meta$type)
heps_from_tumor_counts_boolean=as.data.frame.matrix(heps_from_tumor_counts_boolean)


if(nrow(heps_from_tumor_counts_boolean)==1){
  heps_from_tumor_counts_boolean['1',]=0
  heps_from_tumor_counts_boolean['2',]=0
}else if (nrow(heps_from_tumor_counts_boolean)==2){
  heps_from_tumor_counts_boolean['2',]=0
}



colnames(heps_from_tumor_counts_boolean)=c('Malignant HEPs','Non-Malignant HEPs')


heps_from_tumor_counts_boolean['Level.3',]=colnames(heps_from_tumor_counts_boolean)
heps_from_tumor_counts_boolean['Level.2',]=colnames(heps_from_tumor_counts_boolean)
heps_from_tumor_counts_boolean['Level.1',]='Tumor epithelium'
heps_from_tumor_counts_boolean['Organ',]='Tumor'


#kidney
kidney_boolean=kidney_boolean_whole[row.names(kidney_boolean_whole)%in%pair_of_interest,]
if(is.null(nrow(kidney_boolean))){
  kidney_boolean=rbind(kidney_boolean,0)
}else if(nrow(kidney_boolean)<1){
  kidney_boolean=kidney_boolean_whole[1:2,]
  kidney_boolean[1,]=0
  kidney_boolean[2,]=0
}



kidney_boolean=kidney_boolean[1,]+kidney_boolean[2,]
kidney_boolean=table(number=kidney_boolean,row.names=kidney_meta$cell.type)
kidney_boolean=as.data.frame.matrix(kidney_boolean)
if(nrow(kidney_boolean)==1){
  kidney_boolean['1',]=0
  kidney_boolean['2',]=0
}else if (nrow(kidney_boolean)==2){
  kidney_boolean['2',]=0
}

######################################
kidney_cells=subset(cell_classifier,Organ=='Kidney')

kidney_boolean['Level.3',]=colnames(kidney_boolean)
kidney_boolean['Level.3',]%in%kidney_cells$Level.3
kidney_cells$Level.3%in%kidney_boolean['Level.3',]

kidney_boolean['Level.2',]=kidney_cells[kidney_cells$Level.3==colnames(kidney_boolean),'Level.2']
kidney_boolean['Level.1',]=kidney_cells[kidney_cells$Level.3==colnames(kidney_boolean),'Level.1']
kidney_boolean['Organ',]='Kidney'
####################################


###########################################
#liver

liver_boolean=liver_boolean_whole[row.names(liver_boolean_whole)%in%pair_of_interest,]
if(is.null(nrow(liver_boolean))){
  liver_boolean=rbind(liver_boolean,0)
}else if(nrow(liver_boolean)<1){
  liver_boolean=liver_boolean_whole[1:2,]
  liver_boolean[1,]=0
  liver_boolean[2,]=0
}



liver_boolean=liver_boolean[1,]+liver_boolean[2,]
liver_boolean=table(number=liver_boolean,row.names=liver_meta$cell.type)
liver_boolean=as.data.frame.matrix(liver_boolean)


if(nrow(liver_boolean)==1){
  liver_boolean['1',]=0
  liver_boolean['2',]=0
}else if (nrow(liver_boolean)==2){
  liver_boolean['2',]=0
}


######################################
liver_cells=unique(subset(cell_classifier,Organ=='Liver'))


liver_boolean['Level.3',]=colnames(liver_boolean)
liver_boolean['Level.3',]%in%liver_cells$Level.3
liver_cells$Level.3%in%liver_boolean['Level.3',]

liver_boolean['Level.2',]=liver_cells[liver_cells$Level.3==colnames(liver_boolean),'Level.2']
liver_boolean['Level.1',]=liver_cells[liver_cells$Level.3==colnames(liver_boolean),'Level.1']
liver_boolean['Organ',]='Liver'
####################################
###########################################brain

##########PLUG!!!! PLUG!!!!###################
#brain_boolean=brain_boolean[pair_of_interest,]
brain_boolean=brain_boolean_whole[row.names(brain_boolean_whole)%in%pair_of_interest,]

if(is.null(nrow(brain_boolean))){
  brain_boolean=rbind(brain_boolean,0)
}else if(nrow(brain_boolean)<1){
  brain_boolean=brain_boolean_whole[1:2,]
  brain_boolean[1,]=0
  brain_boolean[2,]=0
}

##############################################
brain_boolean=brain_boolean[1,]+brain_boolean[2,]
brain_boolean=table(number=brain_boolean,row.names=brain_meta$cell.type)
brain_boolean=as.data.frame.matrix(brain_boolean)


if(nrow(brain_boolean)==1){
  brain_boolean['1',]=0
  brain_boolean['2',]=0
}else if (nrow(brain_boolean)==2){
  brain_boolean['2',]=0
}


######################################
brain_cells=unique(subset(cell_classifier,Organ=='Brain'))


brain_boolean['Level.3',]=colnames(brain_boolean)
brain_boolean['Level.3',]%in%brain_cells$Level.3
brain_cells$Level.3%in%brain_boolean['Level.3',]

brain_boolean['Level.2',]=brain_cells[brain_cells$Level.3==colnames(brain_boolean),'Level.2']
brain_boolean['Level.1',]=brain_cells[brain_cells$Level.3==colnames(brain_boolean),'Level.1']
brain_boolean['Organ',]='Brain'
###########################################

#
##########PLUG!!!! PLUG!!!!###################
#brain_boolean=brain_boolean[pair_of_interest,]

lung_boolean=lung_boolean_whole[row.names(lung_boolean_whole)%in%pair_of_interest,]

if(is.null(nrow(lung_boolean))){
  lung_boolean=rbind(lung_boolean,0)
}else if(nrow(lung_boolean)<1){
  lung_boolean=lung_boolean_Whole[1:2,]
  lung_boolean[1,]=0
  lung_boolean[2,]=0
}
##############################################

lung_boolean=lung_boolean[1,]+lung_boolean[2,]
lung_boolean=table(number=lung_boolean,row.names=lung_meta$cell.type)
lung_boolean=as.data.frame.matrix(lung_boolean)


if(nrow(lung_boolean)==1){
  lung_boolean['1',]=0
  lung_boolean['2',]=0
}else if (nrow(lung_boolean)==2){
  lung_boolean['2',]=0
}


######################################
lung_cells=unique(subset(cell_classifier,Organ=='Lung'))


lung_boolean['Level.3',]=colnames(lung_boolean)
lung_boolean['Level.3',]%in%lung_cells$Level.3
lung_cells$Level.3%in%lung_boolean['Level.3',]

lung_boolean['Level.2',]=lung_cells[lung_cells$Level.3==colnames(lung_boolean),'Level.2']
lung_boolean['Level.1',]=lung_cells[lung_cells$Level.3==colnames(lung_boolean),'Level.1']
lung_boolean['Organ',]='Lung'
###########################################
#heart

heart_boolean=heart_boolean_whole[row.names(heart_boolean_whole)%in%pair_of_interest,]

if(is.null(nrow(heart_boolean))){
  heart_boolean=rbind(heart_boolean,0)
}else if(nrow(heart_boolean)<1){
  heart_boolean=heart_booleanheart_boolean_Whole[1:2,]
  heart_boolean[1,]=0
  heart_boolean[2,]=0
}
####

##########PLUG!!!! PLUG!!!!###################
#brain_boolean=brain_boolean[pair_of_interest,]
#heart_boolean=heart_boolean[1:2,]
#heart_boolean[1,]=0
#heart_boolean[2,]=0
##############################################

heart_boolean=heart_boolean[1,]+heart_boolean[2,]
heart_boolean=table(number=heart_boolean,row.names=heart_meta$cell.type)
heart_boolean=as.data.frame.matrix(heart_boolean)


if(nrow(heart_boolean)==1){
  heart_boolean['1',]=0
  heart_boolean['2',]=0
}else if (nrow(heart_boolean)==2){
  heart_boolean['2',]=0
}


######################################
heart_cells=unique(subset(cell_classifier,Organ=='Heart'))


heart_boolean['Level.3',]=colnames(heart_boolean)
heart_boolean['Level.3',]%in%heart_cells$Level.3
heart_cells$Level.3%in%heart_boolean['Level.3',]

heart_boolean['Level.2',]=heart_cells[heart_cells$Level.3==colnames(heart_boolean),'Level.2']
heart_boolean['Level.1',]=heart_cells[heart_cells$Level.3==colnames(heart_boolean),'Level.1']
heart_boolean['Organ',]='Heart'


#colnames(liver_boolean)=paste('Healthy_liver_',colnames(liver_boolean),sep='')
#rm(liver_scrna)
###########################################
#lungs


######overall table
setwd('~/Dropbox/Work/bispecific_markers_project/scrna/')
setwd('D:/Dropbox/Work/bispecific_markers_project/scrna/')
for_heatmap=as.data.frame(t(cbind(heps_from_tumor_counts_boolean,liver_boolean,kidney_boolean,lung_boolean,heart_boolean,brain_boolean)))

write.table(for_heatmap,file=paste0(paste(pair_of_interest,collapse=''),'for_heatmap.txt'))
for_heatmap=read.table(paste0(paste(pair_of_interest,collapse=''),'for_heatmap.txt'))
for_heatmap[for_heatmap$Level.3=='NK cell','Level.3']='NK'
#for_heatmap=as.data.frame(t(for_heatmap))
for_heatmap$sum=apply(for_heatmap[1:3],1,sum)

for_heatmap$X0=for_heatmap$X0/for_heatmap$sum
for_heatmap$X1=for_heatmap$X1/for_heatmap$sum
for_heatmap$X2=for_heatmap$X2/for_heatmap$sum

#for_heatmap=for_heatmap[order(for_heatmap$Level.2,decreasing = T),]

#Level 1
for_heatmap=for_heatmap[order(for_heatmap$X1,for_heatmap$Level.3,decreasing = F),]
#column, Level.1 tissue level aggregation, string column
#Organ. Organ of the origin, string column
#X2= % of the cells expressing both markers, numeric
#X1= % of the cells expressing at least one marker, numeric
level1=ggplot(for_heatmap,aes(y=Level.1,x=Organ))+geom_point(aes(size=X1,fill=X2),pch=21)+
  labs(title=paste0("Co-expression of ",paste(pair_of_interest,collapse=' & ')),
       subtitle='Level 2 tissue aggregation')+
  theme_minimal()+ labs(size = 'At least one \nmarker present',fill='% of cells\nco-expressing\nmarkers')+
  scale_fill_gradient2(mid="white", high="red",limits=c(0,1))+scale_size(range = c(5, 10),limits=c(0,1))+ylab(NULL)

#Level 3
for_heatmap=for_heatmap[order(for_heatmap$X1,for_heatmap$Level.3,decreasing = F),]
level2=ggplot(for_heatmap,aes(y=Level.2,x=Organ))+geom_point(aes(size=X1,fill=X2),pch=21)+
  labs(title=paste0("Co-expression of ",paste(pair_of_interest,collapse=' & ')),
       subtitle='Level 2 tissue aggregation')+
  theme_minimal()+ labs(size = 'At least one \nmarker present',fill='% of cells\nco-expressing\nmarkers')+
  scale_fill_gradient2(mid="white", high="red",limits=c(0,1))+scale_size(range = c(5, 10),limits=c(0,1))+ylab(NULL)

#Level 3
for_heatmap=for_heatmap[order(for_heatmap$X1,for_heatmap$Level.3,decreasing = F),]
level3=ggplot(for_heatmap,aes(y=Level.3,x=Organ))+geom_point(aes(size=X1,fill=X2),pch=21)+
  labs(title=paste0("Co-expression of ",paste(pair_of_interest,collapse=' & ')),
       subtitle='Level 2 tissue aggregation')+
  theme_minimal()+ labs(size = 'At least one \nmarker present',fill='% of cells\nco-expressing\nmarkers')+
  scale_fill_gradient2(mid="white", high="red",limits=c(0,1))+scale_size(range = c(5, 10),limits=c(0,1))+ylab(NULL)

pdf(paste0("Co-expression of ",paste(pair_of_interest,collapse=' & ','.pdf')))
print(level1)
print(level2)
print(level3)
dev.off()
#ggsave(paste0("Co-expression of ",paste(pair_of_interest,collapse=' & ','.pdf')),plot=c(level1,level2,level3))
cat(pair_of_interest,'\n')
}
#################################################