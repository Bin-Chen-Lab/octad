#given a cancer, find its tumor samples, search its normal samples, compute its signature, perform enrichment analysis of its signatures
setwd("~/Documents/stanford/tumor_cell_line/data/pipeline/")
library(DESeq2)
library(tximport)

cancer = "liver hepatocellular carcinoma"
site = "liver"
#dz_phenotype = read.csv("~/Downloads/LIHC_clinicalMatrix", sep = "\t") #could download from treehouse
dz_phenotype = read.csv("~/Downloads/treehouse_public_samples_clinical_metadata.2017-09-11.tsv", sep = "\t")

load("~/Downloads/dz_expr.RData")
GTEX_phenotype =read.csv("~/Downloads/GTEX_phenotype", sep="\t")
#gtex = read.csv("~/Downloads/gtex_RSEM_Hugo_norm_count", sep="\t")
#gtex_subset = gtex[gtex$sample %in% dz_expr$sample, ]

'IDH100 = read.csv("~/Downloads/LGG_IDH100.tsv", sep="\t", stringsAsFactors = F)
IDH200 = read.csv("~/Downloads/LGG_IDH200.tsv", sep="\t", stringsAsFactors = F)
IDH300 = read.csv("~/Downloads/LGG_IDH300.tsv", sep="\t", stringsAsFactors = F)
IDH400 = read.csv("~/Downloads/LGG_IDH400.tsv", sep="\t", stringsAsFactors = F)
IDH1 = rbind(IDH100, IDH200, IDH300, IDH400)
'
dz_expr_patient_reformat = sapply(colnames(dz_expr), function(x){
  paste(unlist(strsplit(x, "\\.")), collapse = "-")
})


dz_expr_patient_tcga_format = sapply(colnames(dz_expr), function(x){
  paste(unlist(strsplit(x, "\\."))[1:3], collapse = "-")
})

sum(dz_expr_patient %in% dz_phenotype$bcr_patient_barcode)

dz_expr_patient_gtex_format = sapply(colnames(dz_expr), function(x){
  paste(unlist(strsplit(x, "\\.")), collapse = "-")
})

sum(dz_expr_patient_gtex_format %in% GTEX_phenotype$Sample)

GTEX_phenotype_ref = GTEX_phenotype$Sample[GTEX_phenotype$body_site_detail..SMTSD. %in% c("Liver")]

dz_tissue = dz_expr[, dz_expr_patient_tcga_format %in%  dz_phenotype$bcr_patient_barcode ] #IDH1$Case.ID
#log2(norm_count+1)
dz_tissue = 2^dz_tissue - 1
ref_tissue = dz_expr[, dz_expr_patient_gtex_format %in% GTEX_phenotype_ref]
ref_tissue = 2^ref_tissue -1

#correlation is stronger while using log
mean(cor(dz_tissue, ref_tissue))

#should we reformat RSEM count
#write.csv(dz_tissue, "~/Downloads/dz_tissue.csv")

#txi_dz = tximport("~/Downloads/dz_tissue.csv", type = "rsem")
#sampleTable <- data.frame(condition = factor(rep(c("A", "B"), each = 3)))
#rownames(sampleTable) <- colnames(txi$counts)
#dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
counts = cbind(dz_tissue, ref_tissue)
rownames(counts) = as.character(dz_expr$sample)
  
coldata = data.frame(sample = colnames(counts) , condition= c(rep("tumor", ncol(dz_tissue)), rep("normal", ncol(ref_tissue))) )

dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = coldata,
                              design= ~ condition)
dds <- DESeq(dds)
rnms <- resultsNames(dds)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition")])
ntd <- normTransform(dds)
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, file="liver/cluster_by_counts.pdf")

rld <- rlog(dds, blind=FALSE)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, file="liver/cluster_by_co_correlation.pdf")

res <- results(dds, contrast=c("condition","tumor","normal"))

plotMA(res, ylim=c(-2,2))
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

res = res[abs(res$log2FoldChange) > 2 & !is.na(res$padj) & res$padj < 0.05, ]
res = res[order(res$log2FoldChange), ]

write.csv(res, paste0("liver/", "liver", "_toil_sig_genes", ".csv")  )

mapping = read.csv("~/Desktop/folder/gene_info_hs.csv")
mapping = mapping[, c("GeneID", "Symbol")]

dz_signature = merge(mapping, data.frame(Symbol = rownames(res), res), by = "Symbol")
dz_signature = dz_signature[order(dz_signature$log2FoldChange), ]
dz_signature = dz_signature[abs(dz_signature$log2FoldChange) > 2 & dz_signature$padj < 0.005, ]
dz_signature$value = dz_signature$log2FoldChange
dz_signature$up_down = ifelse(dz_signature$value > 0, "up", "down")
write.csv(dz_signature, paste0("liver/", "liver", "_toil_sig_genes_mapped", ".csv")  )


