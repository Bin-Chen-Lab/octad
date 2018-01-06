library("GSVA")

load(paste0(dataFolder,"raw/cmpd_sets_", target_type, ".RData"))
cmpdSets = cmpd_sets$cmpd.sets
names(cmpdSets) = cmpd_sets$cmpd.set.names

drug_pred = read.csv(paste0(outputFolder, "/sRGES.csv"))
#drug_pred = read.csv(paste0("~/Documents/stanford/tumor_cell_line/RGES_manuscript/release/data/LIHC/lincs_cancer_sRGES.csv"))

#create a random gene set
random_times = 1000
rgess = matrix(NA, nrow = nrow(drug_pred), ncol = random_times)
for (i in 1:random_times){
  rgess[, i] = sample(drug_pred$sRGES, nrow(rgess))
}
rownames(rgess) = drug_pred$pert_iname
rgess = cbind(rgess, drug_pred$sRGES)

gsea_results = gsva(rgess, cmpdSets, method = "ssgsea",  parallel.sz=8)

gsea_p = apply(gsea_results, 1, function(x){
  sum(x[1:random_times] > x[random_times+1])/random_times
})

gsea_p = data.frame(target = names(gsea_p), p = gsea_p, padj = p.adjust(gsea_p))
gsea_p = gsea_p[order(gsea_p$p), ]

write.csv(gsea_p, paste0(outputFolder, "/enriched_", target_type, ".csv"))


#visualize the top oned
library(limma)
top_target = gsea_p$target[1]
drug_pred$rank = rank(-drug_pred$sRGES)
target_drugs_score = drug_pred$rank[drug_pred$pert_iname %in% cmpdSets[[top_target]]]
pdf(paste0(outputFolder, "/top_enriched_top.pdf"))
  barcodeplot(drug_pred$sRGES, target_drugs_score, main = top_target, xlab = "sRGES")
dev.off()
