#given a cancer, find its tumor samples, search its normal samples, compute its signature, perform enrichment analysis of its signatures
library(DESeq2)
library(edgeR)
#library(tximport)
library(pheatmap)
library(BiocParallel)
library(RColorBrewer)
library(gplots)
#library(Glimma)

####TODO####
#summarise replicates in dz_expr

if (parallel_cores > 1 ){
  register(MulticoreParam(parallel_cores))
}

#load counts and clinical features
####Get Disease Expression Data for the dz parameter####
dz_phenotype <-read.csv(paste0(dataFolder,"raw/treehouse/treehouse_public_samples_clinical_metadata.2017-09-11.tsv"), 
                        sep = "\t", 
                        stringsAsFactors = F)
dz_phenotype$sample_id <- gsub("_","-",dz_phenotype$sample_id)

load(paste0(dataFolder,"raw/treehouse/dz_expr.RData"))
GTEX_phenotype =read.csv(paste0(dataFolder,"raw/treehouse/GTEX_phenotype"), sep="\t", stringsAsFactors = F)
#cancers <-data.frame(table(dz_phenotype$disease)) 

#find cancer sample expression 
#note that dz_expr may have several replicates of samples e.g.
#TCGA-FG-5965-01
#TCGA-FG-5965-02
#TCGA-FG-5965-02.1
#however dz_samples from dz_phenotype only indicate individual samples

dz_samples <-dz_phenotype[dz_phenotype$disease %in% dz,]
dz_tissue <-dz_expr[, colnames(dz_expr) %in%  dz_samples$sample_id ] 
row.names(dz_tissue) <- dz_expr$sample

####TODO summarizing replicate dz_tissue samples####
#1. change colnames by removing rep id e.g.
#TCGA-FG-5965-01
#TCGA-FG-5965-02
#TCGA-FG-5965-02.1
#to TCGA-FG-5965 and average them
  colnames(dz_tissue) <- sapply(colnames(dz_tissue), 
         function(x) paste(unlist(strsplit(x, "-"))[1:3], collapse="-"))
  


# colnames(dz_tissue)<-sapply(colnames(dz_tissue), 
#        function(x) paste(unlist(strsplit(x, "-"))[1:3], collapse="-"))

#some dz_phenotype sample ids are not found in dz_expr samples and vice versa

#annotate dz_samples of tumors with the particular gene mutation
#currently only support single gene mutation 

####mutation tumor vs. non tumor####
if (ncol(dz_tissue) < 4) { stop("few disease tissue samples") }


#outlier detection;
if (remove_outlier == T){
  print("removing outlier dz tissue samples")
  x <-DGEList(counts = round(2^dz_tissue - 1) )
  x <- calcNormFactors(x, method = "TMM")
  lcpm <- cpm(x, log=TRUE)
  pca <- prcomp(t(lcpm))
  #plot(pca$x, pch = 20)
  pc1_z_score = as.numeric(scale(pca$x[,1]))
  pc2_z_score = as.numeric(scale(pca$x[,2]))
  
  outliers = rownames(pca$x)[abs(pc1_z_score) > 3] #
  col.group = rep("black", length(pc1_z_score))
  col.group[rownames(pca$x) %in% outliers] = "red"
  pdf(paste0(outputFolder, "/tissue_mds.pdf"))
    plot(pc1_z_score, pc2_z_score, xlab = "PC1", ylab = "PC2", col = col.group, pch = 20, cex.lab = 1.5, cex.axis = 1.5)
  dev.off()

  #refer to https://www.biostars.org/p/281767/
  dz_tissue = dz_tissue[, !colnames(dz_tissue) %in% outliers]
}

#estimate purity
#may need to use TPM instead of count

if (remove_impure == T){
  print("removing low purity tissue samples")
  sample_purity <-read.csv(paste0(dataFolder,"treehouse_sample_purity.csv"), row.names = 1)
  purity <-sample_purity[colnames(dz_tissue),]
  write.csv(purity, paste0(outputFolder, "/dz_sample_purity.csv"))
  pdf(paste0(outputFolder, "/tumor_purity.pdf"))
    hist(purity, xlab = "purity", cex.lab = 1.5, cex.axis = 1.5, col = "black", main = "")
    abline(v = purity_cutoff, col = "red")
  dev.off()  
  
  dz_tissue <-dz_tissue[, purity > purity_cutoff & !is.na(purity)]
}

####Find and Annotate cancer mutation subtype patient samples####
if (mutation_gene != "" & gdc_project_id != ""){
  library(dplyr)
  print("querying GDC for gene mutation samples")
  gdc_mutated_tumor_pt <-queryGDC(id_mapping_gene_ensembl(mutation_gene), gdc_project_id) %>% filter(gene == 1)
  dz_samples$pt_id <- sapply(dz_samples$sample_id, function(x) paste(unlist(strsplit(x, "-"))[1:3], collapse="-"))
  #summarising dz_samples
  dz_samples <- dz_samples %>% 
    filter(pt_id %in% colnames(dz_tissue)) %>% 
    group_by(pt_id,age_in_years,gender,disease) %>% 
    filter(row_number()==1) %>% ungroup() %>% select(-sample_id)
  dz_samples <- dz_samples %>% mutate(condition = ifelse(dz_samples$pt_id %in% 
                                                           gdc_mutated_tumor_pt$submitter_id, 
                                                         paste0(mutation_gene,"-mutant"),
                                                         paste0("non-",mutation_gene,"-tumor"))) 
  #dz_tissue <- dz_tissue[,colnames(dz_tissue) %in% dz_samples$sample_id]
  
}

#choose top correlated cell lines
#code won't work with these because I stripped the rep number from colnames
#compute_tissue_cell_cor(colnames(dz_tissue))
#compute_tissue_lincs_cell_cor(colnames(dz_tissue))

#find normal sample expression
if (length(site) == 0) {
  print('determining best ref sample tissue site')
  #find all gtex samples
  normal_tissue <-dz_expr[, colnames(dz_expr) %in%  GTEX_phenotype$Sample ] 
  #select varying genes
  iqr_gene <-apply(normal_tissue, 1, IQR)
  varying_genes <-order(iqr_gene, decreasing=T)[1:3000]

  #compare cancer and all normal tissues to infer best reference tissue
  normal_dz_cor <-cor(normal_tissue[varying_genes, ], dz_tissue[varying_genes, ], method = "spearman")
  normal_dz_cor_each <-apply(normal_dz_cor, 1, median)
  
  GTEX_phenotype_cor <-merge(GTEX_phenotype, data.frame(Sample = names(normal_dz_cor_each), cor = as.numeric(normal_dz_cor_each)), by = "Sample")

  write.csv(GTEX_phenotype_cor, paste0(outputFolder, "/GTEX_phenotype_cor.csv"))
  
  #visualize reference tissues
  visualize_top_ref_tissue()
  
  reference_tissue_rank <-aggregate(cor ~ body_site_detail..SMTSD., GTEX_phenotype_cor, median)
  reference_tissue_rank <-reference_tissue_rank[order(reference_tissue_rank$cor, decreasing = T), ]
  
  write.csv(reference_tissue_rank, paste0(outputFolder, "/reference_tissue_rank.csv"))
  site <-reference_tissue_rank$body_site_detail..SMTSD.[1]
  #some normal tissue samples from the best site are not correlated. E.g., in colon ad , Colon - Transverse is the best site, but clearly it has two subgroups
  ref_samples <-intersect(GTEX_phenotype$Sample[tolower(GTEX_phenotype$body_site_detail..SMTSD.) %in% 
                                                  tolower(site)],
                          GTEX_phenotype_cor$Sample[GTEX_phenotype_cor$cor > ref_tissue_cor_cutoff])
  
}else{
  ref_samples <- GTEX_phenotype$Sample[tolower(GTEX_phenotype$body_site_detail..SMTSD.) %in% 
                                         tolower(site)]
}


if (length(ref_samples) < 4) { stop("few reference tissue samples ") }

ref_tissue <- dz_expr[, colnames(dz_expr) %in%  ref_samples ] 

#log2(norm_count+1)
dz_tissue <- 2^dz_tissue - 1
ref_tissue <- 2^ref_tissue -1

#for test(only choose at most 30 samples); it takes time to run DESeq.
#dz_tissue = dz_tissue[, 1:min(30, ncol(dz_tissue))]
#ref_tissue = ref_tissue[, 1:min(30, ncol(ref_tissue))]
#correlation is stronger while using log
#mean(cor(dz_tissue, ref_tissue))

#should we reformat RSEM count
#write.csv(dz_tissue, "~/Downloads/dz_tissue.csv")

#txi_dz = tximport("~/Downloads/dz_tissue.csv", type = "rsem")
#sampleTable <- data.frame(condition = factor(rep(c("A", "B"), each = 3)))
#rownames(sampleTable) <- colnames(txi$counts)
#dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
counts <-cbind(dz_tissue, ref_tissue)
rownames(counts) <-as.character(dz_expr$sample)
coldata <- data.frame(sample = colnames(counts))
coldata <- coldata %>% 
  left_join(dz_samples %>% select(pt_id,condition),
            by=c('sample'='pt_id'))
coldata[is.na(coldata)] <- 'normal'

# coldata <-data.frame(sample = colnames(counts) , 
#                      condition= c(rep("tumor", ncol(dz_tissue)), 
#                                   rep("normal", ncol(ref_tissue))) )

counts <-round(counts)
#hmm.. lots of odd counts?  filtering genes with very large count may miss signficant genes...we may use ruvseqEmpNorm to normalize counts first 
#counts = counts[rowSums(counts) > 0 & rowMax(counts) < 500000, ]
#detect outliers and normalize counts across multiple studies
#need to run the code manually and inspect plots carefully
#need to choose outliers manually
#source("../code/dz/rna_seq_normalization.R") replaced by the function ruvseqEmpNorm
if (normalize_samples == T){
  print('normalizing tissue samples')
  counts <-ruvseqEmpNorm(counts, coldata)
}

####Create Dz and Ref tissue counts DGEList for edgeR or Limma####
x <-DGEList(counts = counts, group = coldata$condition )
x <- calcNormFactors(x, method = "TMM")

cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

#remove lowly expressed genes
keep.exprs <- rowSums(cpm>1) >= min(table(coldata$condition))
x <- x[keep.exprs,, keep.lib.sizes=FALSE]

#unable to run this 
# pdf(paste0(outputFolder, "/tissue_normal_mds.pdf"))
#   #col.group <-coldata$condition
#   #levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
#   #col.group <- as.character(col.group)
#   
#   col.sample <- c("purple","orange")[coldata$condition]
#   
#   plotMDS(x, labels = NULL,  pch = 20, col = col.sample)
#   legend("topleft",fill=c("purple","orange"),legend=levels(coldata$condition))
# dev.off()

#interactive plot
#glMDSPlot(x, labels = colnames(x), groups=x$samples, col = col.group, launch=TRUE)


####DE Method Functions, edgeR, limma, deseq####

if (DE_method == "edgeR"){
  print('computing DE via edgeR')
  dgList <- x
  #order sampleType by normal tissue, tumor w/o specified gene mutation, tumor w/ gene mutation
  #this helps to sort out the design matrix 
  sampleType <-factor(coldata$condition, levels = c('normal', 
                                                       paste0("non-",mutation_gene,"-tumor"),
                                                      paste0(mutation_gene,"-mutant")))
  
  designMat <- model.matrix(~ sampleType)
  
  dgList <- estimateGLMCommonDisp(dgList,design=designMat)
  dgList <- estimateGLMTrendedDisp(dgList,design=designMat)
  dgList <- estimateGLMTagwiseDisp(dgList,design=designMat)
  

  #####calculating DE#####
  #see edgeRUsersGuide section on testing for DE genes for contrast
  fit <- glmFit(dgList,design = designMat)
  
  lrt <- glmLRT(fit,coef=3) #gene mutant tumor vs. ref tissue
  res <- lrt$table
  
  lrt1 <- glmLRT(fit,coef = 2) #non-gene mutant tumor vs. ref tissue
  res1 <- lrt1$table
  
  lrt2 <- glmLRT(fit,contrast=c(0,-1,1)) #gene mutant tumor vs. non-gene mutant tumor
  res2 <- lrt2$table 
  
  res <- lrt$table
  colnames(res) <- c("log2FoldChange", "logCPM", "LR", "pvalue")
  res$padj <- p.adjust(res$pvalue)
  
  res1 <- lrt$table
  colnames(res1) <- c("log2FoldChange", "logCPM", "LR", "pvalue")
  res1$padj <- p.adjust(res1$pvalue)
  
  res2 <- lrt$table
  colnames(res2) <- c("log2FoldChange", "logCPM", "LR", "pvalue")
  res2$padj <- p.adjust(res2$pvalue)
  
  #deciding for significant DE genes
  #deGenes <- decideTestsDGE(lrt,p=0.005)
  #deGenes <- row.names(lrt)[as.logical(deGenes)]
  #plotSmear(lrt,de.tags = deGenes)
  #abline(h=c(-1,1),col=2)
  
}else if (DE_method == "limma"){
  print('computing DE via limma')
  library('Glimma')
  nsamples <- ncol(x)
  col <- brewer.pal(nsamples, "Paired")
  lcpm <- cpm(x, log=TRUE)
  #plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
  #     main="", xlab="")
  #title(main="B. Filtered data", xlab="Log-cpm")
  #abline(v=0, lty=3)
  # for (i in 2:nsamples){
  #   den <- density(lcpm[,i])
  #   lines(den$x, den$y, col=col[i], lwd=2)
  # }
  
  group <-coldata$condition
  col.group <- group
  levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
  col.group <- as.character(col.group)
  #plotMDS(lcpm, labels=coldata$condition, col=col.group)

  ## ----design-----------------------------------------------------------------------------
  design <- model.matrix(~0 + group)
  colnames(design) <- gsub("group", "", colnames(design))
  
  contr.matrix <- makeContrasts(
    TumorvsNon = tumor - normal, 
    levels = colnames(design))
  
  v <- voom(x, design, plot=F)
  vfit <- lmFit(v, design)
  vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
  efit <- eBayes(vfit)
  #plotSA(efit, main="Final model: Mean−variance trend")
  
  tfit <- treat(vfit, lfc=1)
  dt <- decideTests(tfit)
  summary(dt)
  
  tumorvsnormal <- topTreat(tfit, coef=1, n=Inf)
  tumorvsnormal <- tumorvsnormal[order(abs(tumorvsnormal$logFC), decreasing = T),]
  tumorvsnormal.topgenes <- rownames(tumorvsnormal[1:50,])
  'mycol <- colorpanel(1000,"blue","white","red")
  pdf( paste0(outputFolder, "/limma_sig.pdf"))
    heatmap.2(v$E[tumorvsnormal.topgenes,], scale="row",
              labRow=tumorvsnormal.topgenes, labCol=group, 
              col=mycol, trace="none", density.info="none", 
              margin=c(8,6), lhei=c(2,10), dendrogram="column")
  dev.off()'
    
  'df <- as.data.frame(coldata[,c("condition")])
  rownames(df) = (coldata$sample)
  colnames(df) = c("type")
  pheatmap(v$E[tumorvsnormal.topgenes,], cluster_rows=T, show_rownames=T,show_colnames=F, scale="row",col=mycol,
           cluster_cols=T, annotation_col=df, file= paste0(outputFolder, "/limma_sig.pdf"))
  '
  
  res <-tumorvsnormal
  colnames(res) <-c("log2FoldChange", "AveExpr", "t", "pvalue", "padj")
}else{
  print('computing DE via DESeq')
  dds <- DESeqDataSetFromMatrix(countData = round(counts),
                                colData = coldata,
                                design= ~ condition)
  gc()
  if (parallel_cores > 1){
    dds <- DESeq(dds, parallel = T)
  }else{
    dds <- DESeq(dds)
  }
  
  save(dds,file= paste0(outputFolder, "/dds", ".RData"))
  rnms <- resultsNames(dds)
  
  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:50]
  df <- as.data.frame(colData(dds)[,c("condition")])
  rownames(df) <-rownames(colData(dds))
  colnames(df) <-c("type")
  ntd <- normTransform(dds)
  pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
           cluster_cols=FALSE, annotation_col=df, file= paste0(outputFolder, "/deseq_cluster_by_counts.pdf"))
  
  #take time to compute sampel-correlation
  'rld <- rlog(dds, blind=FALSE)
  sampleDists <- dist(t(assay(rld)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors, file= paste0(outputFolder, "/cluster_by_co_correlation.pdf"))
  '
  
  res <- results(dds, contrast=c("condition","tumor","normal"))
}
  
#plotMA(res, ylim=c(-2,2))
#plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

#write.csv(res, paste0(outputFolder, "/dz_sig_genes_all_",DE_method, ".csv")  )

#we missed lots of hits after mapping. need to fix.
#mapping <-read.csv("raw/gene_info_hs.csv")
mapping <- read.csv(paste0(dataFolder,'raw/gene_info_hs.csv'))
mapping <-mapping[, c("GeneID", "Symbol")]

dz_signature <-merge(mapping, data.frame(Symbol = rownames(res), res), by = "Symbol")
dz_signature <-dz_signature[order(dz_signature$log2FoldChange), ]
dz_signature <-dz_signature[abs(dz_signature$log2FoldChange) > dz_fc_threshold & dz_signature$padj < dz_p_threshold & !is.na(dz_signature$Symbol) & !is.na(dz_signature$padj), ]
dz_signature$value <-dz_signature$log2FoldChange
dz_signature$up_down <-ifelse(dz_signature$value > 0, "up", "down")
write.csv(dz_signature, paste0(outputFolder, "/dz_sig_genes_mutant_vs_normal_", DE_method,".csv")  )


dz_signature <-merge(mapping, data.frame(Symbol = rownames(res1), res1), by = "Symbol")
dz_signature <-dz_signature[order(dz_signature$log2FoldChange), ]
dz_signature <-dz_signature[abs(dz_signature$log2FoldChange) > dz_fc_threshold & dz_signature$padj < dz_p_threshold & !is.na(dz_signature$Symbol) & !is.na(dz_signature$padj), ]
dz_signature$value <-dz_signature$log2FoldChange
dz_signature$up_down <-ifelse(dz_signature$value > 0, "up", "down")
write.csv(dz_signature, paste0(outputFolder, "/dz_sig_genes_nonmutantTumor_vs_normal_", DE_method,".csv")  )


dz_signature <-merge(mapping, data.frame(Symbol = rownames(res2), res2), by = "Symbol")
dz_signature <-dz_signature[order(dz_signature$log2FoldChange), ]
dz_signature <-dz_signature[abs(dz_signature$log2FoldChange) > dz_fc_threshold & dz_signature$padj < dz_p_threshold & !is.na(dz_signature$Symbol) & !is.na(dz_signature$padj), ]
dz_signature$value <-dz_signature$log2FoldChange
dz_signature$up_down <-ifelse(dz_signature$value > 0, "up", "down")
write.csv(dz_signature, paste0(outputFolder, "/dz_sig_genes_mutant_vs_nonmutant", DE_method,".csv")  )



#visualize sig genes using lcpm
df <- as.data.frame(coldata[,c("condition")])
rownames(df) <-(coldata$sample)
colnames(df) <-c("type")
mycol <- colorpanel(1000,"blue","white","red")
pheatmap(lcpm[as.character(dz_signature$Symbol),], cluster_rows=T, show_rownames=F,show_colnames=F, scale="row",col=mycol,
         cluster_cols=T, annotation_col=df, file= paste0(outputFolder, "/sig_genes_lcpm.pdf"))

print('Complete dz signature')
#comparing EdgeR, limma, deseq

#In HCC, limma missed the most well-known biomarker "AFP", deseq is super slow, edgeR relatively better in terms of accuracy and speed
