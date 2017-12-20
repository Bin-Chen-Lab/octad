#given a cancer, find its tumor samples, search its normal samples, compute its signature, perform enrichment analysis of its signatures
library(DESeq2)
#library(tximport)
library(pheatmap)
library(BiocParallel)

if (parallel_cores > 1 ){
  register(MulticoreParam(parallel_cores))
}

#dz = "liver hepatocellular carcinoma"
dz_phenotype = read.csv("raw/treehouse/treehouse_public_samples_clinical_metadata.2017-09-11.tsv", sep = "\t", stringsAsFactors = F)
load("raw/treehouse/clinvar.RData")
#reformat clinvar stage data
clinvar$tumor_grade_new = "others"
clinvar$tumor_grade_new[clinvar$tumor_grade %in% c("G2")] = "G2"
clinvar$tumor_grade_new[clinvar$tumor_grade %in% c("G3")] = "G3"

load("raw/treehouse/dz_expr.RData")
GTEX_phenotype =read.csv("raw/treehouse/GTEX_phenotype", sep="\t", stringsAsFactors = F)
cancers = data.frame(table(dz_phenotype$disease))


#dz_expr_patient_reformat = sapply(colnames(dz_expr), function(x){
#  paste(unlist(strsplit(x, "\\.")), collapse = "-")
#})
#colnames(dz_expr) = dz_expr_patient_reformat

#find cancer sample expression
dz_samples = dz_phenotype$sample_id[dz_phenotype$disease %in% dz]

dz = paste0("grade_", dz) 
if (!file.exists(dz)){dir.create(dz)}


dz_tissue = dz_expr[, colnames(dz_expr) %in%  dz_samples ] 

if (ncol(dz_tissue) == 0) { stop() }

#log2(norm_count+1)
dz_tissue = 2^dz_tissue - 1

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
counts = cbind(dz_tissue)
rownames(counts) = as.character(dz_expr$sample)
  
coldata = data.frame(bcr_code = clinvar$bcr_patient_barcode, condition = clinvar$tumor_grade_new) #data.frame(sample = colnames(counts) , condition= c(rep("tumor", ncol(dz_tissue)), rep("normal", ncol(ref_tissue))) )
coldata_1 = merge(coldata, 
                  data.frame(sample = colnames(dz_tissue), 
                             bcr_code = sapply(colnames(dz_tissue), function(x) paste(unlist(strsplit(x, "-"))[1:3], collapse="-"))))
rownames(coldata_1) = coldata_1$sample
counts = counts[, colnames(counts) %in% coldata_1$sample]

if (sum(coldata_1$condition == "G2") < 5 | sum(coldata_1$condition == "G3") < 5){
  next
}

dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = coldata_1[colnames(counts), ],
                              design= ~ condition)

if (parallel_cores > 1){
  dds <- DESeq(dds, parallel = parallel_cores)
}else{
  dds <- DESeq(dds)
}

save(dds,file= paste0(dz, "/dds", ".RData"))

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

res <- results(dds, contrast=c("condition","G3","G2"))

#plotMA(res, ylim=c(-2,2))
#plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

write.csv(res, paste0(dz, "/dz_sig_genes_all", ".csv")  )

#we missed lots of hits after mapping. need to fix.
mapping = read.csv("raw/gene_info_hs.csv")
mapping = mapping[, c("GeneID", "Symbol")]

dz_signature = merge(mapping, data.frame(Symbol = rownames(res), res), by = "Symbol")
dz_signature = dz_signature[order(dz_signature$log2FoldChange), ]
dz_signature = dz_signature[abs(dz_signature$log2FoldChange) > 1.5 & dz_signature$padj < 0.05 & !is.na(dz_signature$padj), ]
dz_signature$value = dz_signature$log2FoldChange
dz_signature$up_down = ifelse(dz_signature$value > 0, "up", "down")
write.csv(dz_signature, paste0(dz, "/dz_sig_genes", ".csv")  )

#valdiate signatures
rnms <- resultsNames(dds)

#select <- order(rowMeans(counts(dds,normalized=TRUE)),
#                decreasing=TRUE)[1:50]
df <- as.data.frame(colData(dds)[,c("condition")])
rownames(df) = rownames(colData(dds))
colnames(df) = c("type")
ntd <- normTransform(dds)
pheatmap(assay(ntd)[dz_signature$Symbol,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, file= paste0(dz, "/cluster_by_counts.pdf"))

