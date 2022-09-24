library(octad)
library(octad.db)
library(dplyr)
library(Seurat)
library(ggrepel)
library(ggplot2)
library(Seurat)
require(gridExtra)
library(octad)
library(octad)
library(EnsDb.Hsapiens.v79)
library(e1071)
library(octad)
library(ROCR)
library(cluster)
library(raster)
library(ggplot2)
library(Seurat)
#Prepare data


setwd('~/Dropbox/Work/bispecific_markers_project/bulk/')
setwd('D:/Dropbox/Work/bispecific_markers_project/bulk/')
surface=read.csv('D:/Dropbox/Work/scRNA/human_compartment_knowledge_full.tsv',sep='\t',header=F)
surface=read.csv('~/Dropbox/Work/scRNA/human_compartment_knowledge_full.tsv',sep='\t',header=F)
table(surface$V4)[order(table(surface$V4),decreasing = T)][1:20]

surface=subset(surface,V4=='Membrane')
surface=surface[grepl('ENSP',surface$V1),]
surface=subset(surface,V7>4)
#surface=subset(surface,V6=='CURATED')
surface_expressed_genes=(unique(surface$V2))
table(surface$V6)




#load results table
result_table=read.table('D:/Dropbox/Work/bispecific_markers_project/bulk/test_SVM_new_4.txt',header=T)
result_table=read.table('~/Dropbox/Work/bispecific_markers_project/bulk/test_SVM_new_4.txt',header=T)

result_table=na.omit(result_table)
#plot(density(result_table$pair_score),lwd=3,col='red',main='Distribution of pair score')
#abline(v=quantile(result_table$pair_score,.90),lwd=4,col='blue')
#result_table=subset(result_table,pair_score>quantile(result_table$pair_score,.90))
#result_table=subset(result_table,pair_score>2)
#result_table=subset(result_table,case_greater=='TRUE_TRUE')
#result_table=result_table[order(result_table$pair_score,decreasing = T),]
#result_table$antigen_1=gsub('-','_',result_table$antigen_1)
#result_table$antigen_2=gsub('-','_',result_table$antigen_2)


head(result_table)
HCC_primary=subset(phenoDF,cancer=='liver hepatocellular carcinoma'&sample.type == 'primary') #select data
case_id=HCC_primary$sample.id #select cases
Healthy=subset(phenoDF,sample.type=='normal'&biopsy.site=='LIVER')
control_id=Healthy$sample.id

cases=loadOctadCounts(case_id,type='tpm',file='D:/Dropbox/binchenlab/work/show/octad.counts.and.tpm.h5')
cases=loadOctadCounts(case_id,type='tpm',file='~/Dropbox/binchenlab/work/show/octad.counts.and.tpm.h5')
cases=as.data.frame(cases)
controls=loadOctadCounts(control_id ,type='tpm',file='D:/Dropbox/binchenlab/work/show/octad.counts.and.tpm.h5')
controls=loadOctadCounts(control_id ,type='tpm',file='~/Dropbox/binchenlab/work/show/octad.counts.and.tpm.h5')
controls=as.data.frame(controls)


#healthy_tissues=subset(phenoDF,sample.type=='normal')
healthy_tissues=subset(phenoDF,sample.type=='normal')  
healthy_tissues=subset(healthy_tissues,grepl('BRAIN',biopsy.site)|biopsy.site=='LIVER'|biopsy.site=='LUNG'|
                         grepl('HEART',biopsy.site)| grepl('KIDNEY',biopsy.site))

healthy_tissues=loadOctadCounts(healthy_tissues$sample.id ,type='tpm',file='D:/Dropbox/binchenlab/work/show/octad.counts.and.tpm.h5')
healthy_tissues=loadOctadCounts(healthy_tissues$sample.id ,type='tpm',file='~/Dropbox/binchenlab/work/show/octad.counts.and.tpm.h5')
healthy_tissues=as.data.frame(healthy_tissues)




#train=cbind(cases,controls,healthy_tissues[sample(ncol(healthy_tissues),ncol(cases))])
train=cbind(cases,controls,healthy_tissues)
#train=cbind(cases,controls)

geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= row.names(train), keytype = "GENEID", columns = c("SYMBOL","GENEID"))

train$GENEID=row.names(train)
train=merge(geneIDs1,train,by='GENEID')
train$GENEID=NULL
train=train[!duplicated(train$SYMBOL),]
row.names(train)=train$SYMBOL
train$SYMBOL=NULL
train=train[surface_expressed_genes,]
#cases=aggregate(cases[-1],by=list(cases$SYMBOL),FUN=mean)
#row.names(cases)=cases$Group.1
#cases$Group.1=NULL
train=na.omit(train)
train$GenestableID=NULL
train$GENEID=NULL

train=as.data.frame(t(train))

train$phenotype=as.factor(c(rep('case',ncol(cases)),rep('control',
                                                        sum(ncol(controls),
                                                            ncol(healthy_tissues)))))

#train$phenotype=as.factor(c(rep('case',ncol(cases)),rep('control',
#                                                        sum(ncol(controls),
#                                                             ncol(healthy_tissues[sample(ncol(healthy_tissues),ncol(cases))])))))



#train$phenotype=as.factor(c(rep('case',ncol(cases)),rep('control',sum(ncol(controls)))))

colnames(train)=gsub('-','_',colnames(train))
#color_vector=c(rep('red',ncol(cases)),rep('blue',ncol(controls)),rep('darkgrey',ncol(healthy_tissues[sample(ncol(healthy_tissues),ncol(cases))])))
color_vector=c(rep('red',ncol(cases)),rep('blue',ncol(controls)),rep('darkgrey',ncol(healthy_tissues)))


#
#PCA
tsne$type <- "others"

HCC_primary=subset(phenoDF,cancer=='liver hepatocellular carcinoma'&sample.type == 'primary') #select data
case_id=HCC_primary$sample.id #select cases
Healthy=subset(phenoDF,sample.type=='normal'&biopsy.site=='LIVER')
control_id=Healthy$sample.id

tsne$type[tsne$sample.id %in% case_id] <- "HCC"
tsne$type[tsne$sample.id %in% control_id] <- "Healthy liver"
tsne$type[tsne$sample.id %in% subset(phenoDF,sample.type=='normal'&grepl('BRAIN',biopsy.site))$sample.id] <- 'Brain'
tsne$type[tsne$sample.id %in% subset(phenoDF,sample.type=='normal'&grepl('LUNG',biopsy.site))$sample.id] <- 'Lung'
tsne$type[tsne$sample.id %in% subset(phenoDF,sample.type=='normal'&grepl('KIDNEY',biopsy.site))$sample.id] <- 'Kidney'
tsne$type[tsne$sample.id %in% subset(phenoDF,sample.type=='normal'&grepl('HEART',biopsy.site))$sample.id] <- 'Heart'

a=subset(tsne,type!='others')
a$X=a$X/5
a$Y=a$Y/5
#plot
PCA <- ggplot(a, aes(X, Y, fill = type)) + geom_point(pch=21,alpha=0.9,size=2)+
  labs(title = paste ('TNSE PLOT'), x= 'TSNE Dim1', y='TSNE Dim2', caption="OCTAD")+
  theme_classic()+labs(fill='Biopsy site')

#ggsave("case_control_map.pdf",height=10,width=10)
#score distribution
result_table=read.table('D:/Dropbox/Work/bispecific_markers_project/bulk/test_SVM_new_4.txt')
result_table=read.table('~/Dropbox/Work/bispecific_markers_project/bulk/test_SVM_new_4.txt')

result_table=na.omit(result_table)

result_table$pair=paste(result_table$antigen_1,result_table$antigen_2,sep=':')
result_table$outlier='NO'
result_table$outlier[result_table$pair_score>quantile(result_table$pair_score,.90)]='YES'
result_table$outlier[result_table$case_greater!='TRUE_TRUE']='NO'
result_table=result_table[order(result_table$pair_score,decreasing = T),]
result_table=subset(result_table,case_greater=='TRUE_TRUE')
# Box plot with dot plot
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)
# Box plot with jittered points
# 0.2 : degree of jitter in x direction
p + geom_jitter(shape=16, position=position_jitter(0.2))

ggplot(data=result_table,aes(x=1,y=pair_score))+geom_boxplot()+ geom_jitter(shape=16, position=position_jitter(0.2))
#fancy volcano plot
ggplot(data=result_table, aes(x=avg_log2FC , y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

#plot(density(result_table$pair_score),lwd=3,col='red',main='Distribution of pair score')
#abline(v=quantile(result_table$pair_score,.90),lwd=4,col='blue')
#GPC3-MUC13
i='CLEC4M'
j='CLEC4G'
plot(type='n',train[,i],train[,j],pch=3,col=color_vector,xlab=i,ylab=j)
points(train[,i][length(color_vector[color_vector!='darkgrey']):ncol(train)],
       train[,j][length(color_vector[color_vector!='darkgrey']):ncol(train)],pch=3,cex=1.3,col='darkgrey')
points(train[,i][1:length(color_vector[color_vector!='darkgrey'])],
       train[,j][1:length(color_vector[color_vector!='darkgrey'])],pch=21,cex=1.3,bg=color_vector)

legend('topleft',legend=c('Tumor','Adjacent','Healthy'),col=c('red','blue','darkgrey'),pch=c(21,21,3),lwd=2,lty=0,pt.bg=c('red','blue','darkgrey'))


#
#
#Bars of sensitivity vs specificity

#Barpots for vital organs for top-10 most frequent genes

#barplots with frequencies

#Volcano plot?


#compose





par(mfrow=c(2,2))

result_table=read.table('D:/Dropbox/Work/bispecific_markers_project/bulk/test_SVM_new_4.txt',header=T)
result_table=read.table('~/Dropbox/Work/bispecific_markers_project/bulk/test_SVM_new_4.txt',header=T)
result_table=result_table[order(result_table$pair_score,decreasing = T),]

result_table = na.omit(result_table)
result_table$significant = "NO"
result_table[result_table$p.adj < 0.01 & result_table$pair_score > 
               quantile(result_table$pair_score,.99) & result_table$case_greater == "TRUE_TRUE", 
             "significant"] = "Enriched in case"
result_table[result_table$p.adj < 0.01 & result_table$pair_score > 
               quantile(result_table$pair_score,.99) & result_table$case_greater == "FALSE_FALSE", 
             "significant"] = "Enriched in control"
result_table$pair_label <- NA
result_table$pair_label[result_table$significant == "Enriched in case"] = apply(result_table[result_table$significant == 
                                                                                               "Enriched in case", ], 1, function(x) paste(x[1:2], 
                                                                                                                                           collapse = ":"))

plot_scores=ggplot(data = result_table, aes(y = pair_score, x = -log10(p.adj), 
                                            col = significant, label = pair_label)) + geom_point() + 
  theme_bw() + ggrepel::geom_text_repel(max.overlaps = 10) + 
  scale_color_manual(values = c("red", "blue", "black"))
##############
#A

#B
plot_scores
#C pca
PCA
#D top1
i='PLVAP'
j='GPC3'

plot(type='n',train[,i],train[,j],pch=3,col=color_vector,xlab=i,ylab=j)
points(train[,i][length(color_vector[color_vector!='darkgrey']):ncol(train)],
       train[,j][length(color_vector[color_vector!='darkgrey']):ncol(train)],pch=3,cex=1.3,col='darkgrey')
points(train[,i][1:length(color_vector[color_vector!='darkgrey'])],
       train[,j][1:length(color_vector[color_vector!='darkgrey'])],pch=21,cex=1.3,bg=color_vector)

legend('topleft',legend=c('Tumor','Adjacent','Healthy'),col=c('red','blue','darkgrey'),pch=c(21,21,3),lwd=2,lty=0,pt.bg=c('red','blue','darkgrey'))


#E muc13-gpc3
i='GPC3'
j='MUC13'
#ggplot(data=train,aes(x=i,y=j))+geom_points()
#g
#ggplot(data=train,aes(x=GPC3,y=MUC13))+geom_point(aes(color=color_vector))
#ggplot(data=train,aes(x=i,y=j))

plot(type='n',train[,i],train[,j],pch=3,col=color_vector,xlab=i,ylab=j)
points(train[,i][length(color_vector[color_vector!='darkgrey']):ncol(train)],
       train[,j][length(color_vector[color_vector!='darkgrey']):ncol(train)],pch=3,cex=1.3,col='darkgrey')
points(train[,i][1:length(color_vector[color_vector!='darkgrey'])],
       train[,j][1:length(color_vector[color_vector!='darkgrey'])],pch=21,cex=1.3,bg=color_vector)

legend('topleft',legend=c('Tumor','Adjacent','Healthy'),col=c('red','blue','darkgrey'),pch=c(21,21,3),lwd=2,lty=0,pt.bg=c('red','blue','darkgrey'))

#F frequency of the markers
#load results table
result_table=read.table('D:/Dropbox/Work/bispecific_markers_project/bulk/test_SVM_new_4.txt',header=T)
result_table=read.table('~/Dropbox/Work/bispecific_markers_project/bulk/test_SVM_new_4.txt',header=T)

result_table=na.omit(result_table)
#plot(density(result_table$pair_score),lwd=3,col='red',main='Distribution of pair score')
#abline(v=quantile(result_table$pair_score,.90),lwd=4,col='blue')
result_table=subset(result_table,pair_score>quantile(result_table$pair_score,.90))
result_table=subset(result_table,pair_score>2)
result_table=subset(result_table,case_greater=='TRUE_TRUE')
result_table=result_table[order(result_table$pair_score,decreasing = T),]
result_table$antigen_1=gsub('-','_',result_table$antigen_1)
result_table$antigen_2=gsub('-','_',result_table$antigen_2)
boxplot=table(c(result_table$antigen_1,result_table$antigen_2))/nrow(result_table)
#boxplot=table(c(result_table$antigen_1,result_table$antigen_2))/nrow(result_table)
boxplot=boxplot[order(boxplot,decreasing = T)][1:30]
barplot=barplot(boxplot,las=2, main='Density of genes across bispecific pairs',cex.names=0.9,col=c('blue','red'),ylim=c(0,0.5))

#D vital organs boxplots
pdf('vital organs boxplots.pdf',width=16,height=8)
for(x in names(boxplot)[1:5]){
  #for_boxplot=list(HCC=train[case_id,x],Liver=healthy_tissues[liver,x],Lung=healthy_tissues[lung,x],
  #                 Kidney=healthy_tissues[kidney,x],Brain=healthy_tissues[brain,x],Heart=healthy_tissues[heart,x],
  #                       Pancreas=healthy_tissues[pancreas,x],Kidney=healthy_tissues[kidney,x],
  #                       Esophagus=healthy_tissues[esophagus,x]
  #                       )
  for_boxplot=list(HCC=train[case_id,x],Liver=healthy_tissues[liver,x],Lung=healthy_tissues[lung,x],
                   Kidney=healthy_tissues[kidney,x],Brain=healthy_tissues[brain,x],Heart=healthy_tissues[heart,x]
  )
  boxplot(for_boxplot,las=2,col=c('red','blue',rep('darkgrey',7)),
          main=paste0('Expression of ',x, ' in different tissues'),
          cex=0.9)
  cat(x,'\n')
}

dev.off()

sapply(LETTERS[1:4], function(x) { 
  plot(rnorm(100))
  fig_label(x, cex=2) 
})
