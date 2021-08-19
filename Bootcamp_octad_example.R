#load the package:
library(octad)


#Example 1. liver hepatocellular carcinoma vs adjacent reference tissues;
#Example 2. breast cancer invasive carcinoma with PIK3 mutation vs reference tissues;
#Example 3. lung adenocarcinoma with amplified MYC gene vs reference tissues;
#Example 4. primary breast cancer invasive carcinoma vs metastatic breast cancer invasive carcinoma unsing only 978 genes expression data for DE;
#Example 5. Compute sRGES score using GEO obtained dataset

#################################
#Example 1. liver hepatocellular carcinoma vs selected reference tissues;
#################################

#1. load phenotype data

head(phenoDF)[1:10] #visualize phenotype table
HCC_primary=subset(phenoDF,cancer=='liver hepatocellular carcinoma'&sample.type == 'primary') #select data
head(HCC_primary)[1:9]
case_id=HCC_primary$sample.id #select cases
write.table(case_id,file='case_id.txt',sep='\t',row.names=F,col.names=F,quote=F)

#2. Compute reference tissue.
#Library pre-requirements and library intallation

#computing reference tissue, by default using small autoencoder, but can use custom expression set,
#by defaulf output=T and outputFolder is empty sending control  corMatrix.csv to working directory
HCC_adjacent=subset(phenoDF,cancer=='liver hepatocellular carcinoma'&sample.type == 'adjacent'&data.source == 'TCGA') #select data

control_id=HCC_adjacent$sample.id #select cases
write.table(control_id,file='control_id.txt',sep='\t',row.names=F,col.names=F,quote=F)


##########OPTIONAL START###############
#visualize relative distance of the computed case and reference ids via precomputed tsne of the OCTAD database
#make a new column called type to state case, control, or others
tsne$type <- "others"
tsne$type[tsne$sample.id %in% case_id] <- "case"
tsne$type[tsne$sample.id %in% control_id] <- "control"

#plot
(p2 <- ggplot(tsne, aes(X, Y, color = type)) + geom_point(alpha = 0.4)+
    labs(title = paste ('TNSE PLOT'), x= 'TSNE Dim1', y='TSNE Dim2', caption="OCTAD")+
    theme_bw())

ggsave("case_control_map.pdf",height=10,width=10)
##########OPTIONAL END###############

#3. Compute DE
#compute differential expression,
#ready expSet='octad.small' or 'octad.whole' using data for LINCS genes or whole octad dataset, IT DOES MATTER!!!, or you can load your own expression data
#by defaulf output=T and outputFolder is empty sending control  corMatrix.csv to working directory
#and you need to source full path to .h5 file containing your TPM expression data. By default files being sourced with package and their names are
# octad.LINCS.h5 for 978 expression signatures and octad.counts.and.tpm.h5
res=diffExp(case_id,control_id,source='octad.whole',output=T,n_topGenes=10000,file='/mnt/research/BigDataBootCamp/day_4_lab/octad/octad.counts.and.tpm.h5')

#Optional in case of absent 'octad.counts.and.tpm.h5'. To download full .h5 file, please reffer to manual.
#res=diffExp(case_id,control_id,source='octad.small',output=T)

#filter result file
res=subset(res,abs(log2FoldChange)>2&padj<0.001)
head(res)

##########OPTIONAL START###############
#Compute  GO enrichment for obtained data
GO = geneEnrich(res$Symbol,database_list=c( "KEGG_2019_Human","GO_Biological_Process_2017","GO_Cellular_Component_2017"),output=T)
head(GO$KEGG_2019_Human)
##########OPTIONAL END###############

#4. Compute sRGES, drug enrichment, cell lines and evaluate results.
#run
sRGES=runsRGES(res,max_gene_size=500,permutations=10000)

#drug enrichment
#sRGES = read.csv('sRGES.csv',stringsAsFactors = F)
octadDrugEnrichment(sRGES = sRGES)

#compute__cell_line. Return most likely lines for case_id
cell_line_computed=computeCellLine(case_id=case_id,returnDF=T)
#process them to evaluate correlation between sRGES and selected cell lines
topLineEval(topline = c('HEPG2'),mysRGES = sRGES)

#################################
#End of example
#################################

#################################
#Example 2. breast cancer invasive carcinoma with PIK3 mutation vs reference tissues;
#################################

#1. load phenotype data
BIC_with_PIK3_primary=subset(phenoDF,cancer=='breast invasive carcinoma'&sample.type == 'primary'&grepl('PIK3CA',mutation_list)) #select data
case_id=BIC_with_PIK3_primary$sample.id #select cases
write.table(case_id,file='case_id.txt',sep='\t',row.names=F,col.names=F,quote=F)
#2. Compute reference tissue.
#computing reference tissue, by default using small autoencoder, but can use custom expression set,
#by defaulf output=T and outputFolder is empty sending control  corMatrix.csv to working directory
control_id=computeRefTissue(case_id,outputFolder='',output=T,adjacent=T,source = "octad",control_size = 50)
write.table(control_id,file='control_id.txt',sep='\t',row.names=F,col.names=F,quote=F)

#3. Compute DE
#compute differential expression,
res=diffExp(case_id,control_id,source='octad.whole',output=T,n_topGenes=10000,file='/mnt/research/BigDataBootCamp/day_4_lab/octad/octad.counts.and.tpm.h5')
#filter result file
res=subset(res,abs(log2FoldChange)>2&padj<0.001)
#4. Compute sRGES, drug enrichment, cell lines and evaluate results.
sRGES=runsRGES(res,max_gene_size=500,permutations=10000)
#drug enrichment
octadDrugEnrichment(sRGES = sRGES)
#compute__cell_line. Return most likely lines for case_id
cell_line_computed=computeCellLine(case_id=case_id,returnDF=T)

#################################
#End of example
#################################

#################################
##Example 3. lung adenocarcinoma with amplified MYC gene vs reference tissues;
#################################

#1. load phenotype data
LA_MYC=subset(phenoDF,cancer=='lung adenocarcinoma'&sample.type == 'primary'&grepl('MYC',gain_list)) #select data
case_id=LA_MYC$sample.id #select cases
write.table(case_id,file='case_id.txt',sep='\t',row.names=F,col.names=F,quote=F)
#2. Compute reference tissue.
#computing reference tissue, by default using small autoencoder, but can use custom expression set,
#by defaulf output=T and outputFolder is empty sending control  corMatrix.csv to working directory
control_id=computeRefTissue(case_id,outputFolder='',output=T,adjacent=T,source = "octad",control_size = 50)
write.table(control_id,file='control_id.txt',sep='\t',row.names=F,col.names=F,quote=F)

#3. Compute DE
#compute differential expression,
res=diffExp(case_id,control_id,source='octad.whole',output=T,n_topGenes=10000,file='/mnt/research/BigDataBootCamp/day_4_lab/octad/octad.counts.and.tpm.h5')
#filter result file
res=subset(res,abs(log2FoldChange)>2&padj<0.001)
#4. Compute sRGES, drug enrichment, cell lines and evaluate results.
sRGES=runsRGES(res,max_gene_size=500,permutations=10000)
#drug enrichment
octadDrugEnrichment(sRGES = sRGES)
#compute__cell_line. Return most likely lines for case_id
cell_line_computed=computeCellLine(case_id=case_id,returnDF=T)

#################################
#End of example
#################################


#################################
#Example 4. primary breast cancer invasive carcinoma vs metastatic breast cancer invasive carcinoma unsing only 978 genes expression data for DE.
#################################

#1. load phenotype data
BIC_primary=subset(phenoDF,cancer=='breast invasive carcinoma'&sample.type == 'primary') #select data
case_id=BIC_primary$sample.id #select cases
write.table(case_id,file='case_id.txt',sep='\t',row.names=F,col.names=F,quote=F)
#2. Compute reference tissue.
#Since we pick metastatic samples, it is available to select them manualy
BIC_metastatic=subset(phenoDF,cancer=='breast invasive carcinoma'&sample.type == 'metastatic') #select data
control_id=BIC_metastatic$sample.id #select cases
write.table(control_id,file='control_id.txt',sep='\t',row.names=F,col.names=F,quote=F)

#3. Compute DE
#compute differential expression, no additional .h5 file is required.
res=diffExp(case_id,control_id,source='octad.small',output=T)

#filter result file
res=subset(res,abs(log2FoldChange)>2&padj<0.001)
#4. Compute sRGES, drug enrichment, cell lines and evaluate results.
sRGES=runsRGES(res,max_gene_size=500,permutations=10000)
#drug enrichment
octadDrugEnrichment(sRGES = sRGES)
#compute__cell_line. Return most likely lines for case_id
cell_line_computed=computeCellLine(case_id=case_id,returnDF=T)

#################################
#End of example
#################################


#################################
#Example 5. Compute sRGES score using GEO obtained dataset
#################################

#In this example we will use HCC vs Adjacent tissue gene counts obtained from GEO via number GSE144269.
download.file(url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144269/suppl/GSE144269_RSEM_GeneCounts.txt.gz',destfile='GSE144269_RSEM_GeneCounts.txt.gz')
data=read.table(gzfile('GSE144269_RSEM_GeneCounts.txt.gz'),header=T,row.names=1)
data=as.matrix(data)
data=log2(data) #log-convert gene expression for edgeR computation
data[is.infinite(data)]=0 #remove infinite numbers

samples=colnames(data) #define the case and control cohorts, A samples in the dataset were obtained from tumors, B samples were obtained from adjacent tissue
case_id=samples[grepl('A_S',samples)]
control_id=samples[grepl('B_S',samples)]

res=diffExp(case_id,control_id,source='side',output=T,n_topGenes=10000,expSet=data,annotate=F) #compute DE
#Since the data already have annotated GENES, there is no use of annotate ENSG data. We can just split the row.names to obtain gene symbols we will use for sRGES computation.
res=cSplit(res,'identifier','|') #Split the identifier
res$Symbol=res$identifier_2 #Set up Symbol column required for sRGES computation
res=subset(res,abs(log2FoldChange)>2&padj<0.01) #Filter the results
sRGES=runsRGES(res,max_gene_size=500,permutations=10000) #conpute sRGES
octadDrugEnrichment(sRGES = sRGES) #compute drug enrichment for obtained scores
#process them to evaluate correlation between sRGES and selected cell lines
rownames(data)<-gsub("\\..*","",rownames(data)) #remove gene alternative transcripts to compute reference cell lines
cell_line_computed=computeCellLine(case_id=case_id,returnDF=T,source='side',expSet=data)
#evaluate top cell line based on the gene expression
topLineEval(topline = c('HEPG2'),mysRGES = sRGES)

#################################
#End of example
#################################
