#dose response/treatment duration analysis

#output_path <- paste(dz, "/all_lincs_score.csv", sep="")
#
output_path = "~/Documents/stanford/breast/release_cmyc/data/oncogene/all_lincs_score_Cyclin_50.csv"
lincs_drug_prediction <- read.csv(output_path)

load("raw/lincs_signatures_cmpd_landmark.RData")

target_expr = lincs_signatures["595", ]

lincs_drug_prediction = merge(lincs_drug_prediction, data.frame(id = names(target_expr), expr = as.numeric(target_expr)), by = "id" )
#should use pert_dose > 0.01
lincs_drug_prediction_subset <- subset(lincs_drug_prediction,  pert_dose > 0 & pert_time %in% c(6, 24))
#pairs that share the same drug and cell id
lincs_drug_prediction_pairs <- merge(lincs_drug_prediction_subset, lincs_drug_prediction_subset, by=c("pert_iname", "cell_id")) 
#x is the reference
lincs_drug_prediction_pairs <- subset(lincs_drug_prediction_pairs, id.x != id.y ) 

length(unique(lincs_drug_prediction_pairs$pert_iname))

dose_independent = lincs_drug_prediction_pairs[lincs_drug_prediction_pairs$pert_dose.x == lincs_drug_prediction_pairs$pert_dose.y,]
dose_independent = dose_independent[dose_independent$pert_time.x > dose_independent$pert_time.y,]

time_independent = lincs_drug_prediction_pairs[lincs_drug_prediction_pairs$pert_time.x == lincs_drug_prediction_pairs$pert_time.y,]
time_independent = time_independent[time_independent$pert_dose.x > time_independent$pert_dose.y, ]


cmpds = unique(time_independent$pert_iname)
expr_norm_counts = NULL
rges_norm_counts = NULL
total_pairs = NULL
p_factor = NULL
p_rges = NULL
for (cmpd in cmpds){
  #print(cmpd)
  time_independent_subset = time_independent[time_independent$pert_iname == as.character(cmpd), ]
  total_pairs = c(total_pairs, nrow(time_independent_subset))
  expr_norm_counts = c(expr_norm_counts, mean((time_independent_subset$expr.x - time_independent_subset$expr.y)/(time_independent_subset$pert_dose.x/time_independent_subset$pert_dose.y))) #/ #sum(time_independent_subset$expr.x < time_independent_subset$expr.y   )/nrow(time_independent_subset))
  rges_norm_counts = c(rges_norm_counts, sum(time_independent_subset$cmap_score.x < time_independent_subset$cmap_score.y   )/nrow(time_independent_subset))
  if (nrow(time_independent_subset) >= 3){
    test = t.test(time_independent_subset$expr.x, time_independent_subset$expr.y, paired = T, alternative = "less")
    p_factor = c(p_factor, test$p.value)
    test = t.test(time_independent_subset$cmap_score.x, time_independent_subset$cmap_score.y, paired = T, alternative = "less")
    p_rges = c(p_rges, test$p.value)
  }else{
    p_factor = c(p_factor, NA)
    p_rges = c(p_rges, NA)
  }
}

results = data.frame(pert_iname = cmpds, expr_norm_counts, rges_norm_counts, total_pairs, p_factor, p_rges)

results = results[order(results$expr_norm_counts, decreasing = T), ]

cmpd = "palbociclib"
time_independent_subset = time_independent[time_independent$pert_iname == as.character(cmpd), ]
t.test(time_independent_subset$expr.x, time_independent_subset$expr.y, paired = T, alternative = "less")
t.test(time_independent_subset$cmap_score.x, time_independent_subset$cmap_score.y, paired = T, alternative = "less")

results[results$pert_iname %in% c("BRD-K73008154"),]

results[results$p_factor < 0.05 & results$p_rges < 0.05 & !is.na(results$p_factor), ]
