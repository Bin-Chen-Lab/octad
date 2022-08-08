library(dplyr)
library(Seurat)
library(ggrepel)
library(ggplot2)
library(Seurat)
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
result_table=read.table('D:/Dropbox/Work/bispecific_markers_project/bulk/test_SVM_new_4.txt')
result_table=read.table('~/Dropbox/Work/bispecific_markers_project/bulk/test_SVM_new_4.txt',header=T)

result_table=na.omit(result_table)
plot(density(result_table$pair_score),lwd=3,col='red',main='Distribution of pair score')
abline(v=quantile(result_table$pair_score,.90),lwd=4,col='blue')
#result_table=subset(result_table,pair_score>quantile(result_table$pair_score,.90))
result_table=subset(result_table,pair_score>2)
result_table=subset(result_table,case_greater=='TRUE_TRUE')
result_table=result_table[order(result_table$pair_score,decreasing = T),]
result_table$antigen_1=gsub('-','_',result_table$antigen_1)
result_table$antigen_2=gsub('-','_',result_table$antigen_2)


head(result_table)


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

geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= row.names(data), keytype = "GENEID", columns = c("SYMBOL","GENEID"))

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

#rm(cases)
#rm(healthy_tissues)
#rm(controls)

#######################################################
#result_table$sens_train=0
#result_table$spec_train=0
#result_table$AUC=0

counter=0
x=1
pdf('result_table_plots_vital_new.pdf')
for (x in 1:nrow(result_table)){
  # for (j in gene_vector){
  i=result_table[x,1]
  j=result_table[x,2]
  counter=counter+1
#  if(i!=j){
    svmfit=svm(as.formula(paste('(phenotype)~',i,'+',j,sep='')),data=train,cost=1,scale=T,cross=5,kernel='radial',probability=TRUE,type='C-classification')
    #svmfit=svm((phenotype)~GPC3+AMBP,data=train,cross=5,type='C-classification',probability=TRUE)
    #fit = svm(factor(y) ~ ., data = dat, scale = FALSE, kernel = "radial", cost = 5)
    #plot(svmfit,train,as.formula(paste(i,'~',j)))
    #plot(svmfit,data=train,GPC3~AMBP)
    #plot(train$GPC3,train$AMBP,col=factor(train$phenotype))
     x.svm.prob <- predict(svmfit, type="prob", newdata=train, probability = TRUE)
     x.svm.prob.rocr <- prediction(as.numeric(as.factor(x.svm.prob)), as.numeric(train$phenotype))
     x.svm.perf <- performance(x.svm.prob.rocr, "tpr","fpr")
    # plot(x.svm.perf, col='red')
    # plot(x.svm.perf, col='red',add=T)
    # x.svm.prob <- predict(svmfit, type="prob", newdata=test, probability = TRUE)
    # x.svm.prob.rocr <- prediction(as.numeric(as.factor(x.svm.prob)), as.numeric(pheno_test))
    # x.svm.perf <- performance(x.svm.prob.rocr, "tpr","fpr")
    # plot(x.svm.perf, col='blue',add=T)
    
    conf=table(predict(svmfit,train),train$phenotype)
    sens_train=conf[1]/(conf[1]+conf[2])
    spec_train=conf[4]/(conf[3]+conf[4])
    result_table[result_table$antigen_1==i&result_table$antigen_2==j,'sens']=sens_train
    result_table[result_table$antigen_1==i&result_table$antigen_2==j,'spec']=spec_train
   # x.svm.perf <- performance(x.svm.prob.rocr, "auc")
    #AUC=x.svm.perf@y.values[[1]]
    # plot(x.svm.perf, col='darkgreen',add=T)
    #result_table[result_table$antigen_1==i&result_table$antigen_2==j,'AUC']=AUC
    
    
    ###markers standalone
    svmfit=svm(as.formula(paste('(phenotype)~',i,sep='')),data=train,cost=1,scale=T,cross=5,kernel='radial',probability=TRUE,type='C-classification')
    conf=table(predict(svmfit,train),train$phenotype)
    sens_train=conf[1]/(conf[1]+conf[2])
    spec_train=conf[4]/(conf[3]+conf[4])
    result_table[result_table$antigen_1==i&result_table$antigen_2==j,'sens_AG1']=sens_train
    result_table[result_table$antigen_1==i&result_table$antigen_2==j,'spec_AG1']=spec_train
    #x.svm.perf <- performance(x.svm.prob.rocr, "auc")
    #AUC=x.svm.perf@y.values[[1]]
    # plot(x.svm.perf, col='darkgreen',add=T)
    #result_table[result_table$antigen_1==i&result_table$antigen_2==j,'AUC_AG1']=AUC
    
    ####################
    svmfit=svm(as.formula(paste('(phenotype)~',j,sep='')),data=train,cost=1,scale=T,cross=5,kernel='radial',probability=TRUE,type='C-classification')
    conf=table(predict(svmfit,train),train$phenotype)
    sens_train=conf[1]/(conf[1]+conf[2])
    spec_train=conf[4]/(conf[3]+conf[4])
    result_table[result_table$antigen_1==i&result_table$antigen_2==j,'sens_AG2']=sens_train
    result_table[result_table$antigen_1==i&result_table$antigen_2==j,'spec_AG2']=spec_train
    #x.svm.perf <- performance(x.svm.prob.rocr, "auc")
    #AUC=x.svm.perf@y.values[[1]]
    # plot(x.svm.perf, col='darkgreen',add=T)
    #result_table[result_table$antigen_1==i&result_table$antigen_2==j,'AUC_AG2']=AUC
    
    ##################################################################################PLOTS
    plot(type='n',train[,i],train[,j],pch=3,col=color_vector,xlab=i,ylab=j)
    points(train[,i][length(color_vector[color_vector!='darkgrey']):ncol(train)],
           train[,j][length(color_vector[color_vector!='darkgrey']):ncol(train)],pch=3,cex=1.3,col='darkgrey')
    points(train[,i][1:length(color_vector[color_vector!='darkgrey'])],
           train[,j][1:length(color_vector[color_vector!='darkgrey'])],pch=21,cex=1.3,bg=color_vector)
    
    legend('topleft',legend=c('Tumor','Adjacent','Healthy'),col=c('red','blue','darkgrey'),pch=c(21,21,3),lwd=0,lty=0,pt.bg=c('red','blue'))
    
    
    
    
    #conf=table(predict(svmfit,test),pheno_test)
    #sens_test=conf[1]/(conf[1]+conf[2])
    #spec_test=conf[4]/(conf[3]+conf[4])
    #result_table[result_table$antigen_1==i&result_table$antigen_2==j,'sens_test']=sens_test
    #result_table[result_table$antigen_1==i&result_table$antigen_2==j,'spec_test']=spec_test
    

    
    
    
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
    # geom_line (aes (y = lm(train[,c(j,i)])))
    #######################################################################
    
    
    
    cat(paste0(round(counter/nrow(result_table)*100,3),'%'),'\n')
    
#}
}
dev.off()
# }
#result_table=na.omit(result_table)


boxplot(c(result_table$spec-result_table$spec_AG1,result_table$spec-result_table$spec_AG2),
        c(result_table$sens-result_table$sens_AG1,result_table$sens-result_table$sens_AG2),
        names=c('Sensitivity','Specificity'),
        main='Change of Sens & Spec for bi- vs monospecific markers')


#write.table(result_table,file='result_table_top_whole_healthy.txt')

result_table=subset(result_table,sens>sens_AG1&sens>sens_AG2&spec>spec_AG1&spec>spec_AG2)


#####################
table(healthy_tissues$biopsy.site)

healthy_tissues=subset(phenoDF,sample.type=='normal')  
#healthy_tissues=subset(healthy_tissues,grepl('COLON',biopsy.site)|biopsy.site=='LIVER'|biopsy.site=='LUNG'|
#                         biopsy.site=='SPLEEN'|biopsy.site=='OVARY'|grepl('PANCREAS',biopsy.site)|
#                         grepl('KIDNEY',biopsy.site)|grepl('ESOPHAGUS',biopsy.site))
healthy_tissues=subset(healthy_tissues,grepl('BRAIN',biopsy.site)|biopsy.site=='LIVER'|biopsy.site=='LUNG'|
                        grepl('HEART',biopsy.site)| grepl('KIDNEY',biopsy.site))

sum(table(healthy_tissues$biopsy.site))


liver=subset(healthy_tissues,biopsy.site=='LIVER')$sample.id
lung=subset(healthy_tissues,biopsy.site=='LUNG')$sample.id
kidney=subset(healthy_tissues,grepl('KIDNEY',biopsy.site))$sample.id
brain=subset(healthy_tissues,grepl('BRAIN',biopsy.site))$sample.id
heart=subset(healthy_tissues,grepl('HEART',biopsy.site))$sample.id




















#spleen=subset(healthy_tissues,biopsy.site=='SPLEEN')$sample.id
#placenta=subset(healthy_tissues,biopsy.site=='OVARY')$sample.id
#pancreas=subset(healthy_tissues,grepl('PANCREAS',biopsy.site))$sample.id
#retina=subset(healthy_tissues,grepl('RETINA',biopsy.site))$sample.id
kidney=subset(healthy_tissues,grepl('KIDNEY',biopsy.site))$sample.id
#esophagus=subset(healthy_tissues,grepl('ESOPHAGUS',biopsy.site))$sample.id
sdf






