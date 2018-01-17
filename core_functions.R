library(estimate)
library(ggplot2)
library("rrcov")
library("RUVSeq")

####TODO####
#need to remove ruvseqEmpNorm should not use normalized counts for the return

####Changes####
#added k param to ruvseq per suggestion of Ke
#gdc_query modified function to allow for querying of tumor project only 
#added compEmpControlGenes function to compute empirical control genes
#added remLowExpr function to remove lowly expressed genes

id_mapping <-function(id, input_format = "hgnc_symbol", output_format = "ensembl_gene_id"){
  #ensembl_gene_id ensembl_transcript_id hgnc_symbol entrezgene
  library(biomaRt)
  
  ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  id_mapped <- getBM(attributes=c(input_format, output_format), filters =
                       input_format, values =id, mart = ensembl)
  return(id_mapped[1,2]) #only return the first hit
}

id_mapping_gene_ensembl <-function(id){
  #ensembl_gene_id ensembl_transcript_id hgnc_symbol entrezgene
  #mapping <- read.csv("~/Documents/GitHub/OCTAD/raw/gencode.v23.annotation.gene.probeMap.csv", stringsAsFactors = F)
  mapping <- read.csv(paste0(dataFolder,'raw/gencode.v23.annotation.gene.probeMap.csv'),stringsAsFactors = F)
  return(mapping$ensembl[mapping$gene == id][1])
}

####GDC API Function####
#Queries GDC website returns the 0 , 1 for genes matching mutation gene or projects
queryGDC <- function(GENE="", PROJECT){
  library(curl)
  
  
  ## examples
  #GENE = id_mapping_gene_ensembl('IDH2')
  #PROJECT = 'TCGA-LGG'
  if (GENE!="") {
    all_gene<-curl_fetch_memory(paste0("https://api.gdc.cancer.gov/analysis/top_mutated_cases_by_gene?fields=submitter_id&pretty=true&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22genes.gene_id%22%2C%22value%22%3A%5B%22", GENE, "%22%5D%7D%7D%5D%7D&format=tsv&size=100000"))
    all_gene.df=read.table(text = rawToChar(all_gene$content), sep = '\t', header = TRUE)    
    all_gene.df$submitter_id=as.character(all_gene.df$submitter_id)    
  }
  
  ## retrieve all IDs with in PROJECT
  all_project<-curl_fetch_memory(paste0("https://api.gdc.cancer.gov/analysis/top_mutated_cases_by_gene?fields=submitter_id&pretty=true&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22",PROJECT,"%22%5D%7D%7D%5D%7D&format=tsv&size=100000"))
  all_project.df=read.table(text = rawToChar(all_project$content), sep = '\t', header = TRUE)    
  all_project.df$submitter_id=as.character(all_project.df$submitter_id)
  
  if(GENE!=""){
  ## combine into grouped matrix
  combined=as.data.frame(unique(append(unique(all_project.df$submitter_id),unique(all_gene.df$submitter_id))))
  colnames(combined)[1]="submitter_id"
  }else{combined=as.data.frame(unique(all_project.df$submitter_id))
        colnames(combined)[1]="submitter_id"}
  
  # add in gene info
  # is 1 
  if(GENE!=''){
  combined$gene=0
  combined[which(combined$submitter_id %in% all_gene.df$submitter_id),"gene"]=1
  }
  
  # add in project info
  combined$project=0
  combined[which(combined$submitter_id %in% all_project.df$submitter_id),"project"]=1  
  
  return(combined)
}



estimatePurity  <- function(expr_matrix){
  #expr_matrix: with gene symbols as row names and samples as colnames
  #return purity score
  samples <- data.frame(NAME = rownames(expr_matrix), Description = NA,  expr_matrix)
  write("", file = "temp_samples.gct")
  write(paste(nrow(samples), "\t", ncol(samples) -2), file = "temp_samples.gct", append = T)
  write("", file = "temp_samples.gct", append = T)
  write.table(samples,  file = "temp_samples.gct", append = T, row.names=F, col.names=T, sep="\t", quote=F)
  
  in.file <- ("temp_samples.gct")
  out.file <- "temp_samples_output.gct"
  estimateScore(in.file, out.file)
  
  estimateScore <- read.delim(out.file, sep = "\t", skip = 3)
  
  return(as.numeric(estimateScore[3, -c(1,2)]))
}

compEmpContGenes <- function(counts, coldata, n_topGenes = 5000){
  set <- newSeqExpressionSet(round(counts),
                             phenoData = data.frame(coldata,row.names= coldata$sample_id))
  
  design <- model.matrix(~ condition, data = pData(set))
  y <- DGEList(counts=counts(set), group =  coldata$condition)
  y <- calcNormFactors(y, method="TMM") #upperquartile generate Inf in the LGG case
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit) #defaults to compare tumor to normal or tumor mutant to normal
  
  top <- topTags(lrt, n=nrow(set))$table
  #n_topGenes <- n_topGenes #5000: assume there are 5000 signficant genes
  
  #based on n_topGenes computing genes with low DE
  #the genes not computed significant is in the empirical set
  empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:n_topGenes]))]
  empirical
}

# ruvseqEmpNorm <- function(counts, coldata, n_topGenes = 5000, k = 1){
#   set <- newSeqExpressionSet(round(counts),
#                              phenoData = data.frame(coldata,row.names= coldata$sample_id))
# 
#   design <- model.matrix(~ condition, data = pData(set))
#   y <- DGEList(counts=counts(set), group =  coldata$condition)
#   y <- calcNormFactors(y, method="TMM") #upperquartile generate Inf in the LGG case
#   y <- estimateGLMCommonDisp(y, design)
#   y <- estimateGLMTagwiseDisp(y, design)
# 
#   fit <- glmFit(y, design)
#   lrt <- glmLRT(fit) #defaults to compare tumor to normal or tumor mutant to normal
# 
#   top <- topTags(lrt, n=nrow(set))$table
#   #n_topGenes <- n_topGenes #5000: assume there are 5000 signficant genes
# 
#   #based on n_topGenes computing genes with low DE
#   #the genes not computed significant is in the empirical set
#   empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:n_topGenes]))]
# 
#   set2 <- RUVg(set, empirical, k=k)
#   #pData(set2)
#   #plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[set2$condition])
#   #plotPCA(set2, col=colors[set2$condition], cex=0.5)
#   colors <- brewer.pal(3, "Set2")
#   
#   pdf(paste0(outputFolder, "/normalized_RNA_Seq_RLE.pdf"))
#     plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[as.factor(coldata$condition)])
#   dev.off()
#   
#   pdf(paste0(outputFolder, "/normalized_RNA_Seq_PCA.pdf"))
#     plotPCA(set2, col=colors[as.factor(coldata$condition)], cex=1.2)
#   dev.off()
#   
#   #output
#   return(set2@assayData$normalizedCounts)
#   
# }

remLowExpr <- function(counts,coldata){
  x <-DGEList(counts = counts, group = coldata$condition )
  #RUVg(x)
  #x <- calcNormFactors(x, method = "TMM")
  #warning output
  #In .calcFactorWeighted(obs = x[, i], ref = x[, refColumn],  ... : NaNs produced
  
  cpm_x <- cpm(x)
  lcpm_x <- cpm(x, log=TRUE)
  
  #remove lowly expressed genes
  ####Possibly where IDH1 may be deleted?####
  keep.exprs <- rowSums(cpm_x>1) >= min(table(coldata$condition))
  keep.exprs
}



compute_tissue_cell_cor <- function(dz_tissue_samples){
  load(paste0(dataFolder,"tissue_cell_line_cor.RData"))
  tumor_cell_cor <- tissue_cell_line_cor[, colnames(tissue_cell_line_cor) %in% dz_tissue_samples]
  #order based on median cor
  tumor_cell_cor_merged <- apply(tumor_cell_cor, 1, median)
  top_cell_lines <- names(head(sort(tumor_cell_cor_merged, decreasing = T), 15))
  tumor_cell_cor <- tumor_cell_cor[top_cell_lines, ]
  #tumor_cell_cor$tumor_type_name = factor(tumor_cell_cor$tumor_type_name, levels = tumor_cell_cor_merged$tumor_type_name)
  
  pdf(paste0(outputFolder, "/top_cell_lines.pdf"))
    par(mar=c(15,4.1,1.1,2.1))
    boxplot(t(tumor_cell_cor), las=2, cex.axis=1)
  dev.off()
}

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


visualize_top_ref_tissue <- function(){
  tissue_ref_cor <- read.csv(paste0(outputFolder, "/GTEX_phenotype_cor.csv"))

  reference_tissue_rank <- aggregate(cor ~ body_site_detail..SMTSD., tissue_ref_cor, median)
  reference_tissue_rank <- reference_tissue_rank[order(reference_tissue_rank$cor, decreasing = T), ]
  
  top_refs <-  reference_tissue_rank[1:10, 1]
  
  tissue_ref_cor <- tissue_ref_cor[tissue_ref_cor$body_site_detail..SMTSD. %in% top_refs, ]
  
  #order based on median cor
  tissue_ref_cor$ref <- factor(tissue_ref_cor$body_site_detail..SMTSD., levels = top_refs)
  
  #margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  pdf(paste0(outputFolder, "/top_reference_tissues.pdf"))
    par(mar=c(12,4.1,4.1,2.1))
    p <- ggplot(tissue_ref_cor, aes(ref, cor))
    print(p +   geom_boxplot(outlier.colour = "grey", notch=F, outlier.shape = NA) + geom_jitter() + theme_bw() + 
            ylab("correlation") +
            xlab("") +  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
                              axis.text.y = element_text(size = 15), axis.title = element_text(size = 20), plot.margin = margin(l=35))
    ) 
  dev.off()  
}


