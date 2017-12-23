######Test Implementation for edgeR#####
library(edgeR)

#reference tutorial
#https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/rnaSeq_DE.pdf

setwd('Documents/GitHub/OCTAD/dz/test code/')

load("ref_tissueSample.RData")
load("dz_tissueSample.RData")
load("edgeR_counts.RData") #only used to merge the gene names


dz_tissue = 2^dz_tissue - 1
ref_tissue = 2^ref_tissue -1
row.names(dz_tissue) = row.names(Counts)
row.names(ref_tissue) = row.names(Counts)
rm(Counts)
counts = cbind(dz_tissue, ref_tissue)



#####normalizing counts RUVSeq methods#####
library("rrcov")
library("RUVSeq")
#library('RColorbrewer')
#library(RColorBrewer)

coldata = data.frame(sample = colnames(counts) , 
                     condition= c(rep("tumor", ncol(dz_tissue)), rep("normal", ncol(ref_tissue))) )
#without normalization

ruvseqEmpNorm <- function(counts,n_topGenes = 5000){
  set <- newSeqExpressionSet(round(counts),
                             phenoData = data.frame(coldata,row.names= coldata$sample))
  
  design <- model.matrix(~ condition, data = pData(set))
  y <- DGEList(counts=counts(set), group =  as.factor(coldata$condition))
  y <- calcNormFactors(y, method="upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef=2)
  
  top <- topTags(lrt, n=nrow(set))$table
  n_topGenes = 5000
  empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:n_topGenes]))]

  set2 <- RUVg(set, empirical, k=1)
  #pData(set2)
  #plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[set2$condition])
  #plotPCA(set2, col=colors[set2$condition], cex=0.5)
  
  #output
  set2@assayData$normalizedCounts
  

}

#library(RColorBrewer)
#colors <- brewer.pal(3, "Set2")
#plotRLE(counts, outline=FALSE, ylim=c(-4, 4), col=colors[set2$condition])

#plotPCA(counts, col=colors[sampleType], cex=0.5)

count_ruvseqNorm <- ruvseqEmpNorm(counts = counts,n_topGenes = 10000)

#plotRLE(count_ruvseqNorm, outline=FALSE, ylim=c(-4, 4), col=colors[set2$condition])
#plotPCA(count_ruvseqNorm, col=colors[sampleType], cex=0.5)
sampleType<- rep("N",ncol(counts))
sampleType[grep("GTEX",colnames(counts))] <- "T"

dgList <- DGEList(counts=count_ruvseqNorm,genes=row.names(counts))

#dgList <- calcNormFactors(dgList,method="TMM")
#plotMDS(dgList)

designMat <- model.matrix(~sampleType)

dgList <- estimateGLMCommonDisp(dgList,design=designMat)
dgList <- estimateGLMTrendedDisp(dgList,design=designMat)
dgList <- estimateGLMTagwiseDisp(dgList,design=designMat)

#plotBCV(dgList) #no idea what I'm looking at...

#####calculating DE#####

fit <- glmFit(dgList,design = designMat)
lrt <- glmLRT(fit)
head(lrt$table)
colnames(lrt$table)
DE_table <- lrt$table
save(DE_table,file="DE_table.RData")

#deciding for significant DE genes
deGenes <- decideTestsDGE(lrt,p=0.005)
deGenes <- row.names(lrt)[as.logical(deGenes)]
plotSmear(lrt,de.tags = deGenes)
abline(h=c(-1,1),col=2)

