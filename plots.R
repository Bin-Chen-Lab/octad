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


boxplots_healthy=head(names(boxplot))


x='GPC3'

##################
#healthy_tissues=subset(phenoDF,sample.type=='normal')
healthy_tissues=subset(phenoDF,sample.type=='normal')  
liver=subset(healthy_tissues,biopsy.site=='LIVER')$sample.id
lung=subset(healthy_tissues,biopsy.site=='LUNG')$sample.id
kidney=subset(healthy_tissues,grepl('KIDNEY',biopsy.site))$sample.id
brain=subset(healthy_tissues,grepl('BRAIN',biopsy.site))$sample.id
heart=subset(healthy_tissues,grepl('HEART',biopsy.site))$sample.id


result_table=subset(result_table,result_table$antigen_1%in%boxplots_healthy|result_table$antigen_2%in%boxplots_healthy)





healthy_tissues=subset(phenoDF,sample.type=='normal')  
healthy_tissues=subset(healthy_tissues,grepl('BRAIN',biopsy.site)|biopsy.site=='LIVER'|biopsy.site=='LUNG'|
                         grepl('HEART',biopsy.site)| grepl('KIDNEY',biopsy.site))
healthy_tissues=loadOctadCounts(healthy_tissues$sample.id ,type='tpm',file='D:/Dropbox/binchenlab/work/show/octad.counts.and.tpm.h5')
healthy_tissues=loadOctadCounts(healthy_tissues$sample.id ,type='tpm',file='~/Dropbox/binchenlab/work/show/octad.counts.and.tpm.h5')
healthy_tissues=as.data.frame(healthy_tissues)



geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= row.names(healthy_tissues), keytype = "GENEID", columns = c("SYMBOL","GENEID"))

healthy_tissues$GENEID=row.names(healthy_tissues)
healthy_tissues=merge(geneIDs1,healthy_tissues,by='GENEID')
healthy_tissues$GENEID=NULL
healthy_tissues=healthy_tissues[!duplicated(healthy_tissues$SYMBOL),]
row.names(healthy_tissues)=healthy_tissues$SYMBOL
healthy_tissues$SYMBOL=NULL
healthy_tissues=healthy_tissues[surface_expressed_genes,]
#cases=aggregate(cases[-1],by=list(cases$SYMBOL),FUN=mean)
#row.names(cases)=cases$Group.1
#cases$Group.1=NULL
healthy_tissues=na.omit(healthy_tissues)
healthy_tissues$GenestableID=NULL
healthy_tissues$GENEID=NULL
healthy_tissues=as.data.frame(t(healthy_tissues))

liver=subset(healthy_tissues,biopsy.site=='LIVER')$sample.id
lung=subset(healthy_tissues,biopsy.site=='LUNG')$sample.id
kidney=subset(healthy_tissues,grepl('KIDNEY',biopsy.site))$sample.id
brain=subset(healthy_tissues,grepl('BRAIN',biopsy.site))$sample.id
heart=subset(healthy_tissues,grepl('HEART',biopsy.site))$sample.id




pdf('per_tissue_boxplots_vital_tissues.pdf',width=14)
for(x in names(boxplot)){
#for_boxplot=list(HCC=train[case_id,x],Liver=healthy_tissues[liver,x],Lung=healthy_tissues[lung,x],
#                 Kidney=healthy_tissues[kidney,x],Brain=healthy_tissues[brain,x],Heart=healthy_tissues[heart,x],
#                       Pancreas=healthy_tissues[pancreas,x],Kidney=healthy_tissues[kidney,x],
#                       Esophagus=healthy_tissues[esophagus,x]
#                       )
for_boxplot=list(HCC=train[case_id,x],Liver=healthy_tissues[liver,x],Lung=healthy_tissues[lung,x],
                 Kidney=healthy_tissues[kidney,x],Brain=healthy_tissues[brain,x],Heart=healthy_tissues[heart,x],
                 )
boxplot(for_boxplot,las=2,col=c('red','blue',rep('darkgrey',7)),
        main=paste0('Expression of ',x, ' in different tissues'),
        cex=0.9)
cat(x,'\n')
}
dev.off()


####################

##################################################################################PLOTS
counter=0
x=1
pdf('result_table_plots_1.pdf')
for (x in 1:nrow(result_table)){
  # for (j in gene_vector){
  i=result_table[x,1]
  j=result_table[x,2]
  counter=counter+1
plot(type='n',train[,i],train[,j],pch=3,col=color_vector,xlab=i,ylab=j)
points(train[,i][length(color_vector[color_vector!='darkgrey']):ncol(train)],
       train[,j][length(color_vector[color_vector!='darkgrey']):ncol(train)],pch=3,cex=1.3,col='darkgrey')
points(train[,i][1:length(color_vector[color_vector!='darkgrey'])],
       train[,j][1:length(color_vector[color_vector!='darkgrey'])],pch=21,cex=1.3,bg=color_vector)

legend('topleft',legend=c('Tumor','Adjacent','Healthy'),col=c('red','blue','darkgrey'),pch=c(21,21,3),lwd=0,lty=0,pt.bg=c('red','blue'))
cat(paste0(round(counter/nrow(result_table)*100,3),'%'),'\n')

#}
}
dev.off()