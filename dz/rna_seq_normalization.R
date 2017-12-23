#detect outliers; visualize PCA and detect manually
#normalize RNA-SEQ data using multiple methods
#Risso D, Ngai J, Speed T and Dudoit S (2014). “Normalization of RNA-seq data using factor analysis of control genes or samples.” Nature Biotechnology, 32(9), pp. 896–902. In press, http://www.nature.com/nbt/journal/v32/n9/full/nbt.2931.html.
#need to check the plots carefully

library("rrcov")
library("RUVSeq")

#some known methods
#looks very sensitive
#pca <- PcaHubert(t(round(counts)),  crit=0.975) 
# length(which(pca@flag == T))
# 
# pca <- prcomp(t(round(counts)))
# plot(pca$sdev)


set <- newSeqExpressionSet(round(counts),
                           phenoData = data.frame(coldata,row.names= coldata$sample))

x = as.factor(coldata$condition)
#without normalization
library(RColorBrewer)
colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)

outliers = "GTEX-T6MN-0011-R4A-SM-5CHSD"
counts = counts[, !colnames(counts) %in% outliers]
coldata = coldata[!coldata$sample %in% outliers,]

set <- newSeqExpressionSet(round(counts),
                           phenoData = data.frame(coldata, row.names= coldata$sample))

plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)

#
'set0 <- betweenLaneNormalization(set, which="upper")
plotRLE(set0, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set0, col=colors[x], cex=1.2)

## ----ruv_spikes, fig.cap="RUVg normalization based on spike-in controls.", fig.subcap=c("RLE plot","PCA plot")----
control.genes <- c("GAPDH", "ACTB")
set1 <- RUVg(set, control.genes, k=1)
pData(set1)
plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set1, col=colors[x], cex=1.2)
'

##emprical
design <- model.matrix(~ condition, data = pData(set))
y <- DGEList(counts=counts(set), group =  as.factor(coldata$condition))
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)

top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]

set2 <- RUVg(set, empirical, k=1)
pData(set2)
pdf(paste0(dz, "/normalized_RNA_Seq_RLE.pdf"))
  plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x])
dev.off()

pdf(paste0(dz, "/normalized_RNA_Seq_RLE.pdf"))
  plotPCA(set2, col=colors[x], cex=1.2)
dev.off()

#output
counts = counts(set2)

