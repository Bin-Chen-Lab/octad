#' @export
####### compute_tissue_lincs_cell_cor #######
compute_tissue_lincs_cell_cor <- function(dz_tissue_samples){
  load(paste0(dataFolder,"tissue_cell_line_cor.RData"))
  ccle_mapping <- read.csv(paste0(dataFolder,"raw/ccle_lincs_mapping.csv"))
  
  tumor_cell_cor <- tissue_cell_line_cor[rownames(tissue_cell_line_cor) %in% ccle_mapping$CCLE.name, colnames(tissue_cell_line_cor) %in% dz_tissue_samples]
  
  #order based on median cor
  tumor_cell_cor_merged <- apply(tumor_cell_cor, 1, median)
  tumor_cell_cor_merged <- merge(data.frame(tumor_cell_cor_merged), ccle_mapping[, c("ccle_cell_line_name", "CCLE.name")], by.x = 0, by.y = "CCLE.name")
  names(tumor_cell_cor_merged) = c("CCLE_name", "cor", "cell_id")
  write.csv(tumor_cell_cor_merged, paste0(outputFolder, "/lincs_cell_lines_cor.csv"))
  
  tumor_cell_cor_merged <- tumor_cell_cor_merged[order(tumor_cell_cor_merged$cor, decreasing = T),]
  top_cell_lines <- tumor_cell_cor_merged$CCLE_name
  tumor_cell_cor <- tumor_cell_cor[top_cell_lines, ]
  #tumor_cell_cor$tumor_type_name = factor(tumor_cell_cor$tumor_type_name, levels = tumor_cell_cor_merged$tumor_type_name)
  
  
  pdf(paste0(outputFolder, "/lincs_cell_lines_cor.pdf"), width = 30)
  par(mar=c(12,4.1,4.1,2.1))
  boxplot(t(tumor_cell_cor), las=2, cex.axis=0.6)
  dev.off()
}