#used to optimize parameters. speed up!
library(DESeq2)
library(edgeR)
#library(tximport)
library(pheatmap)
library(BiocParallel)
library(RColorBrewer)
library(gplots)
#library(Glimma)

if (parallel_cores > 1 ){
  register(MulticoreParam(parallel_cores))
}

if (!file.exists(dz)){dir.create(dz)}


#dz_expr_patient_reformat = sapply(colnames(dz_expr), function(x){
#  paste(unlist(strsplit(x, "\\.")), collapse = "-")
#})
#colnames(dz_expr) = dz_expr_patient_reformat

#find cancer sample expression
dz_samples <-dz_phenotype$sample_id[dz_phenotype$disease %in% dz]
dz_tissue <-dz_expr[, colnames(dz_expr) %in%  dz_samples ] 

#find subtype patient samples
if (mutation_gene != "" & gdc_project_id != ""){
  all_patients <-queryGDC(id_mapping_gene_ensembl(mutation_gene), gdc_project_id)
  dz_tissue_patient_id <-sapply(colnames(dz_tissue), function(x) paste(unlist(strsplit(x, "-"))[1:3], collapse="-"))
  dz_tissue <-dz_tissue[, dz_tissue_patient_id %in% all_patients[all_patients$gene == 1, 1]]
}


if (ncol(dz_tissue) < 4) { stop("few disease tissue samples") }


#outlier detection;
if (remove_outlier == T){
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
  pdf(paste0(dz, "/tissue_mds.pdf"))
    plot(pc1_z_score, pc2_z_score, xlab = "PC1", ylab = "PC2", col = col.group, pch = 20, cex.lab = 1.5, cex.axis = 1.5)
  dev.off()

  #refer to https://www.biostars.org/p/281767/
  dz_tissue = dz_tissue[, !colnames(dz_tissue) %in% outliers]
}

#estimate purity
#may need to use TPM instead of count
if (remove_impure == T){
  sample_purity <-read.csv("treehouse_sample_purity.csv", row.names = 1)
  purity <-sample_purity[colnames(dz_tissue),]
  write.csv(purity, paste0(dz, "/dz_sample_purity.csv"))
  pdf(paste0(dz, "/tumor_purity.pdf"))
    hist(purity, xlab = "purity", cex.lab = 1.5, cex.axis = 1.5, col = "black", main = "")
    abline(v = purity_cutoff, col = "red")
  dev.off()  
  
  dz_tissue <-dz_tissue[, purity > purity_cutoff & !is.na(purity)]
}

#choose top correlated cell lines
compute_tissue_cell_cor(colnames(dz_tissue))
compute_tissue_lincs_cell_cor(colnames(dz_tissue))

#find normal sample expression
if (length(site) == 0) {
  #find all gtex samples
  normal_tissue <-dz_expr[, colnames(dz_expr) %in%  GTEX_phenotype$Sample ] 
  #select varying genes
  iqr_gene <-apply(normal_tissue, 1, IQR)
  varying_genes <-order(iqr_gene, decreasing=T)[1:3000]

  #compare cancer and all normal tissues to infer best reference tissue
  normal_dz_cor <-cor(normal_tissue[varying_genes, ], dz_tissue[varying_genes, ], method = "spearman")
  normal_dz_cor_each <-apply(normal_dz_cor, 1, median)
  
  GTEX_phenotype_cor <-merge(GTEX_phenotype, data.frame(Sample = names(normal_dz_cor_each), cor = as.numeric(normal_dz_cor_each)), by = "Sample")

  write.csv(GTEX_phenotype_cor, paste0(dz, "/GTEX_phenotype_cor.csv"))
  
  #visualize reference tissues
  visualize_top_ref_tissue()
  
  reference_tissue_rank <-aggregate(cor ~ body_site_detail..SMTSD., GTEX_phenotype_cor, median)
  reference_tissue_rank <-reference_tissue_rank[order(reference_tissue_rank$cor, decreasing = T), ]
  
  write.csv(reference_tissue_rank, paste0(dz, "/reference_tissue_rank.csv"))
  site <-reference_tissue_rank$body_site_detail..SMTSD.[1]
}

#some normal tissue samples from the best site are not correlated. E.g., in colon ad , Colon - Transverse is the best site, but clearly it has two subgroups
ref_samples <-intersect(GTEX_phenotype$Sample[tolower(GTEX_phenotype$body_site_detail..SMTSD.) %in% tolower(site)],
                        GTEX_phenotype_cor$Sample[GTEX_phenotype_cor$cor > ref_tissue_cor_cutoff])

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
  
coldata <-data.frame(sample = colnames(counts) , condition= c(rep("tumor", ncol(dz_tissue)), rep("normal", ncol(ref_tissue))) )

counts <-round(counts)
#hmm.. lots of odd counts?  filtering genes with very large count may miss signficant genes...we may use ruvseqEmpNorm to normalize counts first 
#counts = counts[rowSums(counts) > 0 & rowMax(counts) < 500000, ]
#detect outliers and normalize counts across multiple studies
#need to run the code manually and inspect plots carefully
#need to choose outliers manually
#source("../code/dz/rna_seq_normalization.R") replaced by the function ruvseqEmpNorm
if (normalize_samples == T){
  counts <-ruvseqEmpNorm(counts, coldata)
}

x <-DGEList(counts = counts, group = coldata$condition )
x <- calcNormFactors(x, method = "TMM")

cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

#remove lowly expressed genes
keep.exprs <- rowSums(cpm>1) >= min(table(coldata$condition))
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

'pdf(paste0(dz, "/tissue_normal_mds.pdf"))
  #col.group <-coldata$condition
  #levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
  #col.group <- as.character(col.group)
  
  col.sample <- c("purple","orange")[coldata$condition]
  
  plotMDS(x, labels = NULL,  pch = 20, col = col.sample)
  legend("topleft",fill=c("purple","orange"),legend=levels(coldata$condition))
dev.off()
'

#interactive plot
#glMDSPlot(x, labels = colnames(x), groups=x$samples, col = col.group, launch=TRUE)


if (DE_method == "edgeR"){

  dgList <- x
  sampleType <-coldata$condition
  
  designMat <- model.matrix(~ sampleType)
  
  dgList <- estimateGLMCommonDisp(dgList,design=designMat)
  dgList <- estimateGLMTrendedDisp(dgList,design=designMat)
  dgList <- estimateGLMTagwiseDisp(dgList,design=designMat)
  
  #####calculating DE#####
  fit <- glmFit(dgList,design = designMat)
  lrt <- glmLRT(fit)
  #head(lrt$table)
  #colnames(lrt$table)
  
  res <- lrt$table
  colnames(res) <- c("log2FoldChange", "logCPM", "LR", "pvalue")
  res$padj <- p.adjust(res$pvalue)
  
  #deciding for significant DE genes
  #deGenes <- decideTestsDGE(lrt,p=0.005)
  #deGenes <- row.names(lrt)[as.logical(deGenes)]
  #plotSmear(lrt,de.tags = deGenes)
  #abline(h=c(-1,1),col=2)
  
}else if (DE_method == "limma"){

  nsamples <- ncol(x)
  col <- brewer.pal(nsamples, "Paired")
  lcpm <- cpm(x, log=TRUE)
  #plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
  #     main="", xlab="")
  #title(main="B. Filtered data", xlab="Log-cpm")
  #abline(v=0, lty=3)
  #for (i in 2:nsamples){
  #  den <- density(lcpm[,i])
  #  lines(den$x, den$y, col=col[i], lwd=2)
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
  pdf( paste0(dz, "/limma_sig.pdf"))
    heatmap.2(v$E[tumorvsnormal.topgenes,], scale="row",
              labRow=tumorvsnormal.topgenes, labCol=group, 
              col=mycol, trace="none", density.info="none", 
              margin=c(8,6), lhei=c(2,10), dendrogram="column")
  dev.off()'
    
  'df <- as.data.frame(coldata[,c("condition")])
  rownames(df) = (coldata$sample)
  colnames(df) = c("type")
  pheatmap(v$E[tumorvsnormal.topgenes,], cluster_rows=T, show_rownames=T,show_colnames=F, scale="row",col=mycol,
           cluster_cols=T, annotation_col=df, file= paste0(dz, "/limma_sig.pdf"))
  '
  
  res <-tumorvsnormal
  colnames(res) <-c("log2FoldChange", "AveExpr", "t", "pvalue", "padj")
}else{

  dds <- DESeqDataSetFromMatrix(countData = round(counts),
                                colData = coldata,
                                design= ~ condition)
  gc()
  if (parallel_cores > 1){
    dds <- DESeq(dds, parallel = T)
  }else{
    dds <- DESeq(dds)
  }
  
  save(dds,file= paste0(dz, "/dds", ".RData"))
  rnms <- resultsNames(dds)
  
  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:50]
  df <- as.data.frame(colData(dds)[,c("condition")])
  rownames(df) <-rownames(colData(dds))
  colnames(df) <-c("type")
  ntd <- normTransform(dds)
  pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
           cluster_cols=FALSE, annotation_col=df, file= paste0(dz, "/deseq_cluster_by_counts.pdf"))
  
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
           col=colors, file= paste0(dz, "/cluster_by_co_correlation.pdf"))
  '
  
  res <- results(dds, contrast=c("condition","tumor","normal"))
}
  
#plotMA(res, ylim=c(-2,2))
#plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

write.csv(res, paste0(dz, "/dz_sig_genes_all", ".csv")  )

#we missed lots of hits after mapping. need to fix.
mapping <-read.csv("raw/gene_info_hs.csv")
mapping <-mapping[, c("GeneID", "Symbol")]

dz_signature <-merge(mapping, data.frame(Symbol = rownames(res), res), by = "Symbol")
dz_signature <-dz_signature[order(dz_signature$log2FoldChange), ]
dz_signature <-dz_signature[abs(dz_signature$log2FoldChange) > dz_fc_threshold & dz_signature$padj < dz_p_threshold & !is.na(dz_signature$Symbol) & !is.na(dz_signature$padj), ]
dz_signature$value <-dz_signature$log2FoldChange
dz_signature$up_down <-ifelse(dz_signature$value > 0, "up", "down")
write.csv(dz_signature, paste0(dz, "/dz_sig_genes", ".csv")  )

