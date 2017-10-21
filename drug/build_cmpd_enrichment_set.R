#drug enrichment analysis
load("raw/lincs_cmpd_info_mesh_target_sea.RData")

#build enrichment datasets
#mesh
cmpd_sets = list()
target_type = "sea_targets" #sea_targets chembl_targets
lincs_cmpd_info_mesh_target$target = lincs_cmpd_info_mesh_target[, target_type]
lincs_cmpd_info_mesh_target = lincs_cmpd_info_mesh_target[!is.na(lincs_cmpd_info_mesh_target$target), ]

ids = table(unlist(strsplit(paste(lincs_cmpd_info_mesh_target$target, collapse = "\t"), "\t")))
ids = names(ids[ids > 2])

target_cmpd = data.frame()
for(i in 1:nrow(lincs_cmpd_info_mesh_target)){
  target_cmpd = rbind(target_cmpd, data.frame(pert_iname = lincs_cmpd_info_mesh_target$pert_iname[i], target= unlist(strsplit(as.character(lincs_cmpd_info_mesh_target$target[i]), "\t"))))
}
target_cmpd = unique(target_cmpd)

for (i in 1:length(ids)){
  cmpd_sets[["cmpd.sets"]][[i]] = unique(as.vector(target_cmpd$pert_iname[target_cmpd$target == ids[i]]))
  cmpd_sets[["cmpd.set.names"]][[i]] = ids[i]
}

save(cmpd_sets, file = paste0("raw/cmpd_sets_", target_type, ".RData"))




