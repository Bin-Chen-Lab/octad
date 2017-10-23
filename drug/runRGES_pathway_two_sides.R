#an example of running RGES and summarizing RGES across multiple profiles.
#make sure change the workspace and code directory
#setwd("~/Documents/stanford/tumor_cell_line/pipeline/data")

choose_fda_drugs <- FALSE

code_dir <- "../code/drug/"

library("plyr")
library("ggplot2")

fda_drugs = read.csv("raw/repurposing_drugs_20170327.txt", comment.char = "!", header=T, sep="\t")

#load LINCS drug gene expression profiles
load("raw/lincs_signatures_cmpd_landmark_symbol.RData")

#
source(paste(code_dir, "core_functions.R",sep=""))
landmark <- 1
lincs_sig_info <- read.csv("raw/lincs_sig_info.csv")

##############
#compute RGES
#only support landmark genes
gene.list <- rownames(lincs_signatures)


#read pathways
pathway_set_name = "Single_Gene_Perturbations_from_GEO"
pathway_gene_set_up = read.delim(paste0("raw/geneset/", pathway_set_name, "_up.txt"),  sep ="$", header=F, stringsAsFactors = F)
pathway_gene_set_down = read.delim(paste0("raw/geneset/", pathway_set_name, "_down.txt"),  sep ="$", header=F, stringsAsFactors = F)

if (!file.exists(paste0("pathway/", pathway_set_name))) {dir.create(paste0("pathway/", pathway_set_name))}

for (pathway_id in 1:nrow(pathway_gene_set_up)){  #nrow(pathway_gene_set) nrow(pathway_gene_set_up
pathways_up = unlist(strsplit(pathway_gene_set_up[pathway_id,1], "\t\t"))
pathway_name_up = pathways_up[1]

pathways_down = unlist(strsplit(pathway_gene_set_down[pathway_id,1], "\t\t"))
pathway_name_down = pathways_down[1]

pathway_genes_up = as.character(sapply(unlist(strsplit(pathways_up[2], "\t")), function(x) {unlist(strsplit(x, ","))[1]}))
pathway_genes_valid_up = pathway_genes_up[pathway_genes_up %in% gene.list]

pathway_genes_down = as.character(sapply(unlist(strsplit(pathways_down[2], "\t")), function(x) {unlist(strsplit(x, ","))[1]}))
pathway_genes_valid_down = pathway_genes_down[pathway_genes_down %in% gene.list]

if (pathway_name_up != pathway_name_down){
  stop(paste(pathway_id, " name inconsitent"))
}else{
  pathway_name = pathway_name_up
}

print(paste(pathway_id, pathway_name, length(pathway_genes_valid_up), length(pathway_genes_valid_down)))

if (length(pathway_genes_valid_up) <= 4 & length(pathway_genes_valid_down <= 4)){ next}

if (!file.exists(paste0("pathway/", pathway_set_name, "/", pathway_id))) {dir.create(paste0("pathway/", pathway_set_name, "/", pathway_id))}
if (file.exists(paste0("pathway/", pathway_set_name, "/", pathway_id, "/sRGES.csv"))) {next}

output_path <- paste("pathway/", pathway_set_name, "/", pathway_id, "/all_lincs_score.csv", sep="")
sRGES_output_path <- paste("pathway/", pathway_set_name, "/", pathway_id, "/sRGES.csv", sep="")
sRGES_output_path_drug <- paste("pathway/",  pathway_set_name, "/", pathway_id, "/sRGES_drug.csv", sep="")
dz_sig_output_path <- paste("pathway/", pathway_set_name, "/", pathway_id, "/dz_sig_used.csv", sep="")

write(paste(pathway_set_name, pathway_id, pathway_name, length(pathway_genes_valid_up), length(pathway_genes_valid_down) , sep = "\t"), paste0("pathway/", pathway_set_name, "/summary.txt"), append=T)

if (choose_fda_drugs){
  lincs_sig_info <- subset(lincs_sig_info, id %in% colnames(lincs_signatures) & tolower(pert_iname) %in% tolower(fda_drugs$pert_iname))
}else{
  lincs_sig_info <- subset(lincs_sig_info, id %in% colnames(lincs_signatures))
}
#remove duplicate instances
lincs_sig_info <- lincs_sig_info[!duplicated(lincs_sig_info$id),]

sig.ids <- lincs_sig_info$id


##############
#read disease signatures
#dz_signature <- read.csv(paste0(dz, "/dz_sig_genes", ".csv") )
#dz_signature <- read.csv("~/Documents/stanford/breast/release_cmyc/data/myc/drug/dz_signature_lincs.txt", sep = "\t" )

#dz_signature <- dz_signature[dz_signature$GeneID %in% gene.list & abs(dz_signature$value) > 1 & dz_signature$padj < 0.005, ]
#dz_signature$up_down = ifelse(dz_signature$logFC > 1, "up", "down")
dz_signature = rbind(data.frame(GeneID = pathway_genes_valid_up, up_down = "up"), data.frame(GeneID = pathway_genes_valid_down, up_down = "down"))
dz_genes_up <- subset(dz_signature,up_down=="up",select="GeneID")
dz_genes_down <- subset(dz_signature,up_down=="down",select="GeneID")

###############8


#compute RGES
#only choose the top 150 genes
max_gene_size <- 150
if (nrow(dz_genes_up) > max_gene_size){
  dz_genes_up <- data.frame(GeneID= dz_genes_up[1:max_gene_size,])
}
if (nrow(dz_genes_down) > max_gene_size){
  dz_genes_down <- data.frame(GeneID=dz_genes_down[1:max_gene_size,])
}

write.csv(rbind(data.frame(GeneID = dz_genes_up, up_down = "up"), data.frame(GeneID= dz_genes_down, up_down = "down")),  dz_sig_output_path)

dz_cmap_scores <- NULL
count <- 0
for (exp_id in sig.ids) {
  count <- count + 1
 # print(count)
  if (landmark ==1){
    cmap_exp_signature <- data.frame(gene.list,  rank(-1 * lincs_signatures[, as.character(exp_id)], ties.method="random"))    
  }else{
    cmap_exp_signature <- cbind(gene.list,  get.sigs(exp_id))    
  }
  colnames(cmap_exp_signature) <- c("ids","rank")
  dz_cmap_scores <- c(dz_cmap_scores, cmap_score_new(dz_genes_up,dz_genes_down,cmap_exp_signature))
}

#random scores
N_PERMUTATIONS <- 10000 #default 100000
random_sig_ids <- sample(1:ncol(lincs_signatures),N_PERMUTATIONS,replace=T)
count <- 0
random_cmap_scores <- NULL
for (expr_id in random_sig_ids){
  count <- count + 1
  #print(count)
  if (landmark ==1){
    cmap_exp_signature <- data.frame(gene.list,  rank(-1 * lincs_signatures[, as.character(exp_id)], ties.method="random"))    
  }else{
    cmap_exp_signature <- cbind(gene.list,  get.sigs(exp_id))    
  }
  colnames(cmap_exp_signature) <- c("ids","rank")
  
  random_input_signature_genes <- sample(gene.list, (nrow(dz_genes_up)+nrow(dz_genes_down)))
  rand_dz_gene_up <- data.frame(GeneID=random_input_signature_genes[1:nrow(dz_genes_up)])
  rand_dz_gene_down <- data.frame(GeneID=random_input_signature_genes[(nrow(dz_genes_up)+1):length(random_input_signature_genes)])
  random_cmap_scores <- c(random_cmap_scores, cmap_score_new(rand_dz_gene_up,rand_dz_gene_down,cmap_exp_signature))
}

p <- sapply(dz_cmap_scores, function(score){
  sum(random_cmap_scores < score)/length(random_cmap_scores)
})

padj <- p.adjust(p, "fdr")
results <- data.frame(id = sig.ids, cmap_score = dz_cmap_scores, p, padj)

results <- merge(results, lincs_sig_info, by = "id")
results <- results[order(results$cmap_score),]
write.csv(results, output_path)

####################
#summarize RGES
lincs_drug_prediction <- read.csv(output_path)

#should use pert_dose > 0.01
lincs_drug_prediction_subset <- subset(lincs_drug_prediction,  pert_dose > 0 & pert_time %in% c(6, 24))
#pairs that share the same drug and cell id
lincs_drug_prediction_pairs <- merge(lincs_drug_prediction_subset, lincs_drug_prediction_subset, by=c("pert_iname", "cell_id")) 
#x is the reference
lincs_drug_prediction_pairs <- subset(lincs_drug_prediction_pairs, id.x != id.y & pert_time.x == 24 & pert_dose.x == 10) #, select <- c("cmap_score.x", "cmap_score.y", "pert_dose.y", "pert_time.y"))

#difference of RGES to the reference 
lincs_drug_prediction_pairs$cmap_diff <- lincs_drug_prediction_pairs$cmap_score.x - lincs_drug_prediction_pairs$cmap_score.y
lincs_drug_prediction_pairs$dose <- round(log(lincs_drug_prediction_pairs$pert_dose.y, 2), 1)

#estimate difference
lincs_drug_prediction_pairs$dose_bin <- ifelse(lincs_drug_prediction_pairs$pert_dose.y < 10, "low", "high")
diff <- tapply(lincs_drug_prediction_pairs$cmap_diff, paste(lincs_drug_prediction_pairs$dose_bin, lincs_drug_prediction_pairs$pert_time.y), mean)

#ignore weighting cell lines
#lincs_cell_line_weight <- read.csv("lincs_cell_line_weight.csv")
#pred <- merge(lincs_drug_prediction, lincs_cell_line_weight, by.x="cell_id", by.y="lincs_cell_id")
pred <- lincs_drug_prediction
pred$cor <- 1
pred$RGES <- sapply(1:nrow(pred), function(id){getsRGES(pred$cmap_score[id], pred$cor[id], pred$pert_dose[id], pred$pert_time[id], diff, max(pred$cor))})

cmpd_freq <- table(pred$pert_iname)
pred <- subset(pred, pert_iname %in% names(cmpd_freq[cmpd_freq > 0]))

pred_merged <- ddply(pred,  .(pert_iname),  summarise,
                     mean = mean(RGES),
                     n = length(RGES),
                     median = median(RGES),
                     sd = sd(RGES))
pred_merged$sRGES <- pred_merged$mean
pred_merged <- pred_merged[order(pred_merged$sRGES), ]
write.csv(pred_merged, sRGES_output_path)

#limit to FDA approved drugs
if (choose_fda_drugs){
  pred_merged_drug <- merge(fda_drugs, pred_merged, by = "pert_iname")
  pred_merged_drug <- pred_merged_drug[order(pred_merged_drug$sRGES), ]
  write.csv(pred_merged_drug, sRGES_output_path_drug)
}

}