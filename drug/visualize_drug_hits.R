#visualize reversal
library(pheatmap)
library("gplots")

#CMAP score output
output_path <- paste(dz, "/all_lincs_score.csv", sep="")
lincs_drug_prediction = read.csv(output_path)

dz_gene_used_output_path <- paste(dz, "/dz_sig_used.csv", sep="")
dz_gene_ids = read.csv(dz_gene_used_output_path)

dz_signature <- read.csv(paste0(dz, "/dz_sig_genes", ".csv") )
dz_signature = dz_signature[dz_signature$GeneID %in% dz_gene_ids$GeneID, ]
dz_sig = subset(dz_signature, select=c("GeneID", "log2FoldChange"))

#drugs to visualize
sRGES = read.csv(paste(dz, "/sRGES.csv", sep=""))
#choose top common drugs
top_drugs = as.character(unique(sRGES$pert_iname[sRGES$n > 1 & !sRGES$pert_iname %in% sRGES$pert_iname[grep("BRD-|SA-", sRGES$pert_iname)]][1:20]))

##########
#visualized gene reversed
#only pick the signatures from close to the median
#cell_lines = read.csv(paste("../table/", cancer, "_cell_lines.csv", sep=""))
#lincs_drug_prediction_subset = subset(lincs_drug_prediction, cell_id %in% cell_lines$LINCS, select=c("id", "RGES", "pert_iname", "pert_dose", "pert_time"))
lincs_drug_prediction$RGES = lincs_drug_prediction$cmap_score
lincs_drug_prediction_subset = lincs_drug_prediction[lincs_drug_prediction$pert_iname %in% top_drugs,]
###selecting median still sounds weird... let's keep all signatures 
drug_cmap_score = aggregate(RGES ~ pert_iname, lincs_drug_prediction_subset, median)
drug_instances_median  = merge(lincs_drug_prediction_subset, drug_cmap_score, by = c("pert_iname"))
drug_instances_median$diff = abs(drug_instances_median$RGES.x - drug_instances_median$RGES.y) #cmap_score.y is the median
drug_instances_min_diff = aggregate(diff ~ pert_iname, drug_instances_median, min)
drug_instances_select = merge(drug_instances_median, drug_instances_min_diff, by=c("pert_iname", "diff"))
drug_instances_select = drug_instances_select[!duplicated(drug_instances_select$pert_iname), ]
sig_id_selects = drug_instances_select$id

#sig_id_selects = lincs_drug_prediction_subset$id

if (landmark == 1){
  load("raw/lincs_signatures_cmpd_landmark.RData")
}else{
  load("raw/lincs_signatures_cmpd_landmark_GSE92742.RData")
}

drug_dz_signature = merge(dz_sig, data.frame(GeneID = rownames(lincs_signatures), lincs_signatures[, as.character(sig_id_selects)]),  by="GeneID", suffixes='')


#########################
###
#visualize the reversed gene expression
#reorder drug dz signatures
gene_ids = drug_dz_signature$GeneID
drug_dz_signature_rank = drug_dz_signature[,-1]
for (i in 1:ncol(drug_dz_signature_rank)){
  drug_dz_signature_rank[,i] = rank(-1 * drug_dz_signature_rank[,i] ) #highly expressed genes ranked on the top
}
gene_ids_rank <- gene_ids[order(drug_dz_signature_rank[,1])]
drug_dz_signature_rank <- drug_dz_signature_rank[order(drug_dz_signature_rank[,1]),] #order by disease expression

col_sorted = sort(cor(drug_dz_signature_rank, method="spearman")["log2FoldChange",-1])    
drug_dz_signature_rank = drug_dz_signature_rank[,c("log2FoldChange", names(col_sorted))]

drug_names = sapply(2:ncol(drug_dz_signature_rank), function(id){
  lincs_drug_prediction$pert_iname[paste("X",lincs_drug_prediction$id, sep="") == names(drug_dz_signature_rank)[id]][1]
})

pdf(paste(dz, "/lincs_reverse_expression.pdf", sep=""))
#colPal <- bluered(100)
colPal = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256))
par(mar=c(13, 6, 2, 0.5))
axiscolor = sapply(c(dz, as.character(drug_names)), function(name){
  if (name == dz){
    "black"
  }else if (name %in% ""){
    "black"
  }else{
    "black"
  }
})
image(t(drug_dz_signature_rank), col=colPal,   axes=F, srt=45)
axis(1,  at= seq(0,1,length.out=ncol( drug_dz_signature_rank ) ), labels= FALSE)
text(x = seq(0,1,length.out=ncol( drug_dz_signature_rank ) ), c(-0.05),
     labels = c( dz,as.character(drug_names)),col=axiscolor, srt = 45, pos=2,offset=0.05, xpd = TRUE, cex=0.6)
dev.off()
