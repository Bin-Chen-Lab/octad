library("GSVA")

load(paste0("raw/cmpd_sets_", target_type, ".RData"))
cmpdSets = cmpd_sets$cmpd.sets
names(cmpdSets) = cmpd_sets$cmpd.set.names

drug_pred = read.csv(paste0(dz, "/sRGES.csv"))
#drug_pred = read.csv(paste0("~/Documents/stanford/tumor_cell_line/RGES_manuscript/release/data/LIHC/lincs_cancer_sRGES.csv"))

#create a random set
random_times = 1000
rgess = matrix(NA, nrow = nrow(drug_pred), ncol = random_times)
for (i in 1:random_times){
  rgess[, i] = sample(drug_pred$sRGES, nrow(rgess))
}
rownames(rgess) = drug_pred$pert_iname
rgess = cbind(rgess, drug_pred$sRGES)

gsea_results = gsva(rgess, cmpdSets, method = "ssgsea",  parallel.sz=0)

gsea_p = apply(gsea_results, 1, function(x){
  sum(x[1:random_times] > x[random_times+1])/random_times
})

gsea_p = data.frame(target = names(gsea_p), p = gsea_p, padj = p.adjust(gsea_p))
gsea_p = gsea_p[order(gsea_p$p), ]

write.csv(gsea_p, paste0(dz, "/enriched_", target_type, ".csv"))
