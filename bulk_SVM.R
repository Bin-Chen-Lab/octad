
library(octad)
library(EnsDb.Hsapiens.v79)
library(e1071)
library(octad)
library(ROCR)
library(cluster)
library(raster)
library(ggplot2)
library(Seurat)

setwd('~/Dropbox/Work/SVM')
#surface=read.csv('~/Dropbox/Work/cell.surface.marker.csv',header=T,sep=';')
#
#

surface=read.csv('D:/Dropbox/Work/scRNA/human_compartment_knowledge_full.tsv',sep='\t',header=F)
surface=read.csv('~/Dropbox/Work/scRNA/human_compartment_knowledge_full.tsv',sep='\t',header=F)
table(surface$V4)[order(table(surface$V4),decreasing = T)][1:20]

surface=subset(surface,V4=='Membrane')
surface=surface[grepl('ENSP',surface$V1),]
surface=subset(surface,V7>4)
#surface=subset(surface,V6=='CURATED')
surface_expressed_genes=(unique(surface$V2))
table(surface$V6)


#load expression data for raw counts or tpm values.
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


# 1. Convert from ensembl.gene to gene.symbol
#ensembl.genes <- genes
#cases$GenestableID=row.names(cases)
#cases=merge(surface[1:2],cases,by='GenestableID')
#cases$GenestableID=NULL
#cases=aggregate(cases[-1],by=list(cases$Gene.name),FUN=mean)
#row.names(cases)=cases$Group.1
#cases$Group.1=NULL


#controls$GenestableID=row.names(controls)
#controls=merge(surface[1:2],controls,by='GenestableID')
#controls$GenestableID=NULL
#controls=aggregate(controls[-1],by=list(controls$Gene.name),FUN=mean)
#row.names(controls)=controls$Group.1
#controls$Group.1=NULL



geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= row.names(cases), keytype = "GENEID", columns = c("SYMBOL","GENEID"))


cases$GENEID=row.names(cases)
cases=merge(geneIDs1,cases,by='GENEID')
cases$GENEID=NULL
cases=cases[!duplicated(cases$SYMBOL),]
row.names(cases)=cases$SYMBOL
cases$SYMBOL=NULL
cases=cases[surface_expressed_genes,]
#cases=aggregate(cases[-1],by=list(cases$SYMBOL),FUN=mean)
#row.names(cases)=cases$Group.1
#cases$Group.1=NULL
cases=na.omit(cases)
cases$GenestableID=NULL

controls$GENEID=row.names(controls)
controls=merge(geneIDs1,controls,by='GENEID')
controls$GENEID=NULL
controls=controls[!duplicated(controls$SYMBOL),]
row.names(controls)=controls$SYMBOL
controls$SYMBOL=NULL
controls=controls[surface_expressed_genes,]
controls=na.omit(controls)
controls$GenestableID=NULL


#controls=aggregate(controls[-1],by=list(controls$SYMBOL),FUN=mean)
#row.names(controls)=controls$Group.1
#controls$Group.1=NULL
annotation=data.frame(sample=c(colnames(cases),colnames(controls)),phenotype=c(rep('cancer',length(colnames(cases))),rep('control',length(colnames(controls)))))
annotation$phenotype=as.factor(annotation$phenotype)


expression=DGEList(counts=as.matrix(cbind(cases,controls)),group=annotation$phenotype)
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
DE$table$padj=p.adjust(DE$table$PValue,method='BH')
DE=subset(DE$table,padj<0.01&abs(logFC)>1.3)
#DE=subset(DE,padj<0.01&abs(logFC)>1.3)

################################ prepare data for SVM
cases=cases[surface$V2,]
cases=cases[row.names(DE),]
#cases['phenotype',]='case'
#cases['phenotype',]=1
controls=controls[surface$V2,]
controls=controls[row.names(DE),]
#controls['phenotype',]='control'
#controls['phenotype',]=2

train_case_vector=sample(1:length(colnames(cases)),length(colnames(cases))*0.75,replace=F)
train_control_vector=sample(1:length(colnames(controls)),length(colnames(controls))*0.75,replace=F)

cases[is.na(cases)]=0
controls[is.na(controls)]=0


train_cases=cases[train_case_vector]
train_control=controls[train_control_vector]
test_cases=cases[-train_case_vector]
test_control=controls[-train_control_vector]

train=as.data.frame(t(cbind(train_cases,train_control)))
#train
test=as.data.frame(t(cbind(test_cases,test_control)))

pheno_train=as.factor(c(rep('case',ncol(train_cases)),rep('control',ncol(train_control))))
pheno_test=as.factor(c(rep('case',ncol(test_cases)),rep('control',ncol(test_control))))
#train$phenotype=as.numeric(train$phenotype)
#test$phenotype=as.numeric(test$phenotype)
#train$phenotype=as.factor(train$phenotype)
#test$phenotype=as.factor(test$phenotype)

#all_pairs=t(combn(colnames(train),2))

#result_table=data.frame(antigen_1=all_pairs[,1],
#                        antigen_2=all_pairs[,2],
#                        sens_train=0,spec_train=0,
#                        sens_test=0,spec_test=0,
#                        AUC=0,distance=0,spread=0,
#                        pair_score=0)

i=colnames(train)[1]
j=colnames(train)[2]
counter=0

gene_vector=c(colnames(train))
all_pairs=t(combn(gene_vector,2))
result_table=data.frame(antigen_1=all_pairs[,1],
                        antigen_2=all_pairs[,2],
                        distance=0,spread=0,angle_cos=0,
                        pair_score=0,case_greater=FALSE)
a=Sys.time()
counter=0
#pdf('svmfit.pdf')
for (x in 1:nrow(result_table)){
 # for (j in gene_vector){
  i=result_table[x,1]
  j=result_table[x,2]
    counter=counter+1
    if(i!=j){
#svmfit=svm(as.formula(paste('pheno_train~',i,'+',j,sep='')),data=train,kernel='radial',cross=5,type='C-classification',probability=TRUE)
#svmfit=svm(as.formula(paste('factor(phenotype)~',i,'+',j,sep='')),
#          data=train,cost=1,scale=T,cross=5,kernel='radial',probability=TRUE,type='C-classification')
#svmfit=svm((phenotype)~GPC3+AMBP,data=train,cross=5,type='C-classification',probability=TRUE)
#fit = svm(factor(y) ~ ., data = dat, scale = FALSE, kernel = "radial", cost = 5)


##########define centers and dist

distance=sqrt(sum((pam(train[which(pheno_train=='case'),c(j,i)], 1)$medoids-
                       pam(train[which(pheno_train=='control'),c(j,i)], 1)$medoids)^2)) #distance between clusters
result_table[result_table$antigen_1==i&result_table$antigen_2==j,'distance']=distance

#silhouette

temp_dist=dist(train[,c(j,i)])
silh=(mean(silhouette(as.numeric(pheno_train),temp_dist)[,3]))+1
result_table[result_table$antigen_1==i&result_table$antigen_2==j,'spread']=silh

#angle coefficient
table_temp=rbind(pam(train[which(pheno_train=='case'),c(j,i)], 1)$medoids,
                 pam(train[which(pheno_train=='control'),c(j,i)], 1)$medoids)
#Boolean logic to check the position of clusters
  result_table[result_table$antigen_1==i&result_table$antigen_2==j,'case_greater']=paste0(table_temp[1,]>table_temp[2,],collapse ='_')
  
m1=lm(table_temp[,i]~table_temp[,j])$coefficients[2]
cos_to_45=1/(1+(abs((1-m1)/(1+m1*1)))^2)
result_table[result_table$antigen_1==i&result_table$antigen_2==j,'angle_cos']=cos_to_45
#score
pair_score=cos_to_45*silh*distance
result_table[result_table$antigen_1==i&result_table$antigen_2==j,'pair_score']=pair_score



#plot(silh,col=c('red','blue'))

#temp_dist=dist(x)
#silh=mean(silhouette(as.numeric(pheno_train[1:500]),temp_dist)[,3])
#plot(silh,col=c('red','blue'))
#plot(train[,i]~train[,j],col=pheno_train)
#train[which(pheno_train=='case'),c(i,j)]
#train[1:10,c(i,j)]
#train[pheno_train[pheno_train=='control'],c(i,j)]
#plot(train[,i]~train[,j],col=pheno_train)
#points(pam(train[which(pheno_train=='case'),c(j,i)], 1)$medoids,pch=16,cex=4,col='red')
#points(pam(train[which(pheno_train=='control'),c(j,i)], 1)$medoids,pch=16,cex=4,col='blue')

#lines(1:20,1:20)
#abline(coef=c(4,m1))

#silhouette(pheno_train,dist(train[,c(j,i)]))


#library(cluster)
#df <- data.frame(X = rnorm(100, 0), Y = rpois(100, 2))
#plot(df$X, df$Y)
#points(pam(df, 1)$medoids, pch = 16, col = "red")
#
#ggplot(data=train[,c(j,i)],aes(AMBP,GPC3))+geom_point()+geom_smooth (method = "loess") +
#  geom_line (aes (y = lm(train[,c(j,i)])))
#######################################################################


#plot(svmfit,train,as.formula(paste(i,'~',j)))
#plot(svmfit,data=train,GPC3~AMBP)
#plot(train$GPC3,train$AMBP,col=factor(train$phenotype))

#x.svm.prob <- predict(svmfit, type="prob", newdata=train, probability = TRUE)
#x.svm.prob.rocr <- prediction(as.numeric(as.factor(x.svm.prob)), as.numeric(pheno_train))
#x.svm.perf <- performance(x.svm.prob.rocr, "tpr","fpr")
#plot(x.svm.perf, col='red')
#plot(x.svm.perf, col='red',add=T)
#x.svm.prob <- predict(svmfit, type="prob", newdata=test, probability = TRUE)
#x.svm.prob.rocr <- prediction(as.numeric(as.factor(x.svm.prob)), as.numeric(pheno_test))
#x.svm.perf <- performance(x.svm.prob.rocr, "tpr","fpr")
#plot(x.svm.perf, col='blue',add=T)

#conf=table(predict(svmfit,train),pheno_train)
#sens_train=conf[1]/(conf[1]+conf[2])
#spec_train=conf[4]/(conf[3]+conf[4])
#result_table[result_table$antigen_1==i&result_table$antigen_2==j,'sens_train']=sens_train
#result_table[result_table$antigen_1==i&result_table$antigen_2==j,'spec_train']=spec_train

#conf=table(predict(svmfit,test),pheno_test)
#sens_test=conf[1]/(conf[1]+conf[2])
#spec_test=conf[4]/(conf[3]+conf[4])
#result_table[result_table$antigen_1==i&result_table$antigen_2==j,'sens_test']=sens_test
#result_table[result_table$antigen_1==i&result_table$antigen_2==j,'spec_test']=spec_test

#x.svm.perf <- performance(x.svm.prob.rocr, "auc")
#AUC=x.svm.perf@y.values[[1]]
#plot(x.svm.perf, col='darkgreen',add=T)
#result_table[result_table$antigen_1==i&result_table$antigen_2==j,'AUC']=AUC

cat(paste0(round(counter/nrow(result_table)*100,3),'%'),'\n')

    }}
#  }
a-Sys.time()
result_table=na.omit(result_table)


dev.off()
#write.table(result_table,'test_SVM_new_4.txt')
result_table=read.table('test_SVM_new_4.txt')

result_table=na.omit(result_table)
result_table=subset(result_table,case_greater=='TRUE_TRUE')
result_table=subset(result_table,pair_score>quantile(result_table$pair_score,.90))

plot(density(result_table$pair_score))

result_table=na.omit(result_table)

#result_table=subset(result_table,angle_cos>0.7)
plot(density(result_table$pair_score))


result_table$antigen_1=gsub('-','_',result_table$antigen_1)
result_table$antigen_2=gsub('-','_',result_table$antigen_2)
#colnames(train)=gsub('-','_',colnames(train))

pdf('test_result_5.pdf')
for(i in 1:nrow(result_table)){
  svmfit=svm(as.formula(paste('pheno_train~',result_table[i,1],'+',result_table[i,2],sep='')),data=train,kernel='radial',cross=5,type='C-classification',probability=TRUE)
  plot(svmfit,train,as.formula(paste(result_table[i,2],'~',result_table[i,1])))
  cat(i/nrow(result_table)*100,'\n')
}
dev.off()


svmfit=svm(as.formula(paste('pheno_train~',i,'+',j,sep='')),data=train,kernel='radial',cross=5,type='C-classification',probability=TRUE)
plot(svmfit,train,as.formula(paste(j,'~',i)))
###################################################
#scRNA
###################################################
load('~/Dropbox/Work/scRNA/HCC_final_pipeline_scRNA/tumor_HEPS.Rds')
load('~/Dropbox/Work/scRNA/HCC_final_pipeline_scRNA/tumor_T_cells_counts.Rds')
load('~/Dropbox/Work/scRNA/HCC_final_pipeline_scRNA/healthy_HEPS.Rds')
load('~/Dropbox/Work/scRNA/HCC_final_pipeline_scRNA/healthy_T_cells_counts.Rds')


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
heps <- FindClusters(heps, resolution = 0.5)
head(Idents(heps), 5)
heps <- RunUMAP(heps, dims = 1:30)

rm(healthy_heps_counts)
rm(tumor_heps_counts)

FeaturePlot(object = heps, features = c("MUC13", "GPC3"), cols= c("grey", "red", "blue"),blend = TRUE)
tumor_markers <- FindMarkers(heps, ident.1 = levels(heps)[grepl('T',levels(heps))], 
                             ident.2 = levels(heps)[!grepl('T',levels(heps))], only.pos = TRUE)

#tumor_markers=subset(tumor_markers,pct.1-pct.2>0.5)
tumor_markers_surface=tumor_markers[row.names(tumor_markers)%in%unique(c(result_table$antigen_1,result_table$antigen_2)),]
# view results
head(monocyte.de.markers)

subset(result_table,result_table$antigen_1%in%row.names(tumor_markers_surface)|result_table$antigen_2%in%row.names(tumor_markers_surface))


pdf('cd4+gpc3+ambp.pdf',height=12,width=12)
CombinePlots(plots = list(
  #  FeaturePlot(heps, features = c('CD4'),label=T),
  #  FeaturePlot(heps, features = c('AMBP'),label=T),
  FeaturePlot(heps, features = c('GPC3'),label=T),
  DimPlot(heps, reduction = "umap",label=T,group.by='pheno')
))
dev.off()

