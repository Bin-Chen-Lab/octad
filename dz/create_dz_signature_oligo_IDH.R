#given a cancer, find its tumor samples, search its normal samples, compute its signature, perform enrichment analysis of its signatures
#oligo: only select IDH mutated vs IDH wild type

library(DESeq2)
library(tximport)
library(pheatmap)
library(edgeR)
dz = "glioma"
site = NULL

mapping = read.csv("raw/gene_info_hs.csv")
mapping = mapping[, c("GeneID", "Symbol")]

dz_phenotype = read.csv("raw/treehouse/treehouse_public_samples_clinical_metadata.2017-09-11.tsv", sep = "\t", stringsAsFactors = F)
load("raw/treehouse/dz_expr.RData")
GTEX_phenotype =read.csv("raw/treehouse/GTEX_phenotype", sep="\t", stringsAsFactors = F)
cancers = data.frame(table(dz_phenotype$disease))

if (!file.exists(dz)){dir.create(dz)}

IDH100 = read.csv("~/Downloads/LGG_IDH100.tsv", sep="\t", stringsAsFactors = F)
IDH200 = read.csv("~/Downloads/LGG_IDH200.tsv", sep="\t", stringsAsFactors = F)
IDH300 = read.csv("~/Downloads/LGG_IDH300.tsv", sep="\t", stringsAsFactors = F)
IDH400 = read.csv("~/Downloads/LGG_IDH400.tsv", sep="\t", stringsAsFactors = F)
IDH1 = rbind(IDH100, IDH200, IDH300, IDH400)

dz_expr_patient_reformat = sapply(colnames(dz_expr), function(x){
  paste(unlist(strsplit(x, "\\.")), collapse = "-")
})

colnames(dz_expr) = dz_expr_patient_reformat

dz_expr_patient_tcga_format = sapply(colnames(dz_expr), function(x){
  paste(unlist(strsplit(x, "-"))[1:3], collapse = "-")
})

sum(dz_expr_patient_tcga_format %in% IDH1$Case.ID)

dz_expr_patient_gtex_format = sapply(colnames(dz_expr), function(x){
  paste(unlist(strsplit(x, "\\.")), collapse = "-")
})

dz_phenotype$bcr_patient = sapply(dz_phenotype$sample_id, function(x){
  paste(unlist(strsplit(x, "-"))[1:3], collapse = "-")
})


#find cancer sample expression
dz_samples = dz_phenotype$sample_id[ dz_phenotype$disease %in% dz]
dz_samples_IDH_mutant = dz_phenotype$sample_id[dz_phenotype$bcr_patient %in% IDH1$Case.ID]
dz_samples_IDH_wild = dz_phenotype$sample_id[!dz_phenotype$bcr_patient %in% IDH1$Case.ID & dz_phenotype$disease %in% dz]

dz_tissue = dz_expr[, colnames(dz_expr) %in%  dz_samples ] 

#find normal sample expression
if (length(site) == 0) {
  #find all gtex samples
  normal_tissue = dz_expr[, colnames(dz_expr) %in%  GTEX_phenotype$Sample ] 
  #select varying genes
  iqr_gene = apply(normal_tissue, 1, IQR)
  varying_genes = order(iqr_gene, decreasing=T)[1:3000]

  #compare cancer and all normal tissues to infer best reference tissue
  normal_dz_cor = cor(normal_tissue[varying_genes, ], dz_tissue[varying_genes, ], method = "spearman")
  normal_dz_cor_each = apply(normal_dz_cor, 1, median)
  
  GTEX_phenotype_cor = merge(GTEX_phenotype, data.frame(Sample = names(normal_dz_cor_each), cor = as.numeric(normal_dz_cor_each)), by = "Sample")
  
  reference_tissue_rank = aggregate(cor ~ body_site_detail..SMTSD., GTEX_phenotype_cor, median)
  reference_tissue_rank = reference_tissue_rank[order(reference_tissue_rank$cor, decreasing = T), ]
  
  write.csv(reference_tissue_rank, paste0(dz, "/reference_tissue_rank.csv"))
  site = reference_tissue_rank$body_site_detail..SMTSD.[1]
}

site = "Brain - Cortex" #suggested from Anders

ref_samples = GTEX_phenotype$Sample[tolower(GTEX_phenotype$body_site_detail..SMTSD.) %in% tolower(site)]
ref_tissue = dz_expr[, colnames(dz_expr) %in%  ref_samples ] 

#log2(norm_count+1)
#hmm, the total number of counts in TCGA and GTEX are different
dz_tissue = 2^dz_tissue - 1
ref_tissue = 2^ref_tissue -1

hist(apply(dz_tissue, 2, sum))
summary(apply(dz_tissue, 2, sum))
hist(apply(ref_tissue, 2, sum))
summary(apply(ref_tissue, 2, sum))
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
counts = cbind(dz_tissue, ref_tissue)
rownames(counts) = as.character(dz_expr$sample)

#only keep protein-code genes... lots of non-coding here

counts = counts[rownames(counts) %in% mapping$Symbol, ]
dz_conditions = sapply(colnames(dz_tissue), function(x){
  if (x %in% dz_samples_IDH_mutant){
    "IDH_mutant"
  }else if (x %in% dz_samples_IDH_wild){
    "IDH_wild"
  }
})
coldata = data.frame(sample = colnames(counts) , condition= c(dz_conditions,  rep("normal", ncol(ref_tissue))) )

dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = coldata,
                              design= ~ condition)

dds <- DESeq(dds, parallel=TRUE)



rnms <- resultsNames(dds)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]
df <- as.data.frame(colData(dds)[,c("condition")])
rownames(df) = rownames(colData(dds))
colnames(df) = c("type")
ntd <- normTransform(dds)
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, file= paste0(dz, "/cluster_by_counts.pdf"))

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
#we missed lots of hits after mapping. need to fix.

save(dds, file = paste0(dz, "/dds.RData"))
load(paste0(dz, "/dds.RData"))
#MUTATION
res <- results(dds, contrast=c("condition","IDH_mutant","normal"))
write.csv(res, paste0(dz, "/dz_sig_genes_all_IDH_mutant_normal", ".csv")  )

dz_signature = merge(mapping, data.frame(Symbol = rownames(res), res), by = "Symbol")
dz_signature = dz_signature[order(dz_signature$log2FoldChange), ]
dz_signature = dz_signature[abs(dz_signature$log2FoldChange) > 2 & dz_signature$padj < 0.005, ]
dz_signature$value = dz_signature$log2FoldChange
dz_signature$up_down = ifelse(dz_signature$value > 0, "up", "down")
write.csv(dz_signature, paste0(dz, "/dz_sig_genes_IDH_mutant_normal", ".csv")  )

#wild
res <- results(dds, contrast=c("condition","IDH_wild", "normal"))
write.csv(res, paste0(dz, "/dz_sig_genes_all_IDH_wild_normal", ".csv")  )

dz_signature = merge(mapping, data.frame(Symbol = rownames(res), res), by = "Symbol")
dz_signature = dz_signature[order(dz_signature$log2FoldChange), ]
dz_signature = dz_signature[abs(dz_signature$log2FoldChange) > 2 & dz_signature$padj < 0.005, ]
dz_signature$value = dz_signature$log2FoldChange
dz_signature$up_down = ifelse(dz_signature$value > 0, "up", "down")
write.csv(dz_signature, paste0(dz, "/dz_sig_genes_IDH_wild_normal", ".csv")  )


#wild vs mutation
res <- results(dds, contrast=c("condition","IDH_mutant","IDH_wild"))
write.csv(res, paste0(dz, "/dz_sig_genes_all_IDH_mutant_wild", ".csv")  )

dz_signature = merge(mapping, data.frame(Symbol = rownames(res), res), by = "Symbol")
dz_signature = dz_signature[order(dz_signature$log2FoldChange), ]
dz_signature = dz_signature[abs(dz_signature$log2FoldChange) > 2 & dz_signature$padj < 0.005, ]
dz_signature$value = dz_signature$log2FoldChange
dz_signature$up_down = ifelse(dz_signature$value > 0, "up", "down")
write.csv(dz_signature, paste0(dz, "/dz_sig_genes_IDH_mutant_wild", ".csv")  )

