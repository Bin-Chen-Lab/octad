#given a cancer, find its tumor samples, search its normal samples, compute its signature, perform enrichment analysis of its signatures
#library(DESeq2) #moved this to DESeq2 method
library(edgeR)
#library(tximport)
library(pheatmap)
library(BiocParallel)
library(RColorBrewer)
library(gplots)
library(dplyr)
#library(Glimma)
test = 0
####DataFrames####
#dz_expr : DF of all log2 expression counts of dz and ref tissues
  #dz_tissue : subset of dz_expr containing the cancer tissue of interest
    # this df may remove computed impure or outliers 
#dz_phenotype : feature data of dz in dz_expr
  #dz_tissue_phenotype : feature data of dz in selected dz_tissues 
    # this df will not remove computed impure or outliers but will annotate them as such
#GTEX_phenotype : feature data of normal tissue samples in dz_expr
  #ref_samples_phenotype : subset of ref samples chosen to make comparisons with dz tissues
#counts : absolute value (2 ^ expr - 1) of dz_tissue and dz_phenotype
#counts_phenotype : phenotype data of counts
  #counts$condition : this should be the condition you want to run your contrast in e.g.
    #tumor vs. normal
    #mutant gene tumor vs. normal, mutant gene tumor vs. non-mutant gene tumor, non-mutant gene tumor vs. normal


####Changes####
#all_pheno is the dataframe for phenodata of dz_expr 
#GTEX_phenotype and dz_phenotype no longer needed as all_pheno is more comprehensive
#remove outliers made into a function as part of core_functions

####TODO####
#make function for GTEx tissue selection, may need to embed w/ visualize function
  #parameters
    #random tissues
    #varying genes
    #correlation cutoffs
    #visualize function 
#update methods for LIMMA and DESeq


if (parallel_cores > 1 ){
  register(MulticoreParam(parallel_cores))
}

#load counts and clinical features
#turn test on if you're testing code the test dz_expr Sample set has ~10000 samples 
#test = 0
if(test == 1)
{
  load(paste0(dataFolder,"raw/treehouse/dz_exprSamp.RData"))
}else{load(paste0(dataFolder,"raw/treehouse/dz_exprUncomp.RData"))}

all_pheno <- read.csv(paste0(dataFolder,"/raw/treehouse/TcgaTargetGTEX_phenotype.txt"),sep="\t")
all_pheno <- all_pheno %>% filter(X_primary_site != '')

# dz_phenotype <-read.csv(paste0(dataFolder,"raw/treehouse/treehouse_public_samples_clinical_metadata.2017-09-11.tsv"), 
#                         sep = "\t", 
#                         stringsAsFactors = F)
# dz_phenotype$sample_id <- gsub("_","-",dz_phenotype$sample_id)


#For TCGA The fourth values after dash should be samples
#note that dz_expr may have several replicates of samples e.g.
#TCGA-FG-5965-01
#TCGA-FG-5965-02
#TCGA-FG-5965-02.1

GTEX_phenotype =read.csv(paste0(dataFolder,"raw/treehouse/GTEX_phenotype"), sep="\t", stringsAsFactors = F)
#cancers <-data.frame(table(dz_phenotype$disease)) 
all_pheno$patient_id <- sapply(all_pheno$sample, function(x) paste(unlist(strsplit(as.character(x), "-"))[1:3], collapse="-"))
rownames(dz_expr) <- dz_expr$sample
dz_expr <- dz_expr[,colnames(dz_expr) %in% all_pheno$sample]



####Select Tumor####
#there are multiple sample types not only tumor,
#check sample type column to make sure you're getting what you need

dz_tissue_phenotype <- all_pheno %>% filter(toupper(primary.disease.or.tissue) == toupper(dz),
                                            grepl('tumor',tolower(X_sample_type)))
dz_tissue_phenotype$sample_id <- dz_tissue_phenotype$sample


dz_tissue <- dz_expr[, colnames(dz_expr) %in% dz_tissue_phenotype$sample ] 
#row.names(dz_tissue) <- dz_expr$sample
dz_tissue_phenotype <- dz_tissue_phenotype %>% 
  filter(sample %in% colnames(dz_tissue)) #maintain consistency
#some dz_phenotype sample ids are not found in dz_expr samples and vice versa

#annotate dz_samples of tumors with the particular gene mutation
#currently only support single gene mutation 

####Data Cleaning : Compute and Remove Outlier, Impure (optional)####
if (ncol(dz_tissue) < 4) { stop("few disease tissue samples") }


#outlier detection moved to core_functions 

if (remove_outlier == T){
  print("removing outlier dz tissue samples")
  outliers <- detectOutlier(dz_tissue,3)
  numOutliers = length(outliers)
  ifelse(numOutliers == 0,
         print('no outliers found'),
         print(paste0(numOutliers," outliers removed")))
  dz_tissue = dz_tissue[, !colnames(dz_tissue) %in% outliers]
  rm(numOutliers,outliers)
  dz_tissue_phenotype <- dz_tissue_phenotype %>% 
    filter(sample %in% colnames(dz_tissue)) #maintain consistency
}
  
#estimate purity
#may need to use TPM instead of count

if (remove_impure == T){
  print("removing low purity tissue samples")
  sample_purity <-read.csv(paste0(dataFolder,"treehouse_sample_purity.csv"), row.names = 1)
  purity <-sample_purity[colnames(dz_tissue),]
  write.csv(purity, paste0(outputFolder, "/dz_sample_purity.csv"))
  pdf(paste0(outputFolder, "/tumor_purity.pdf"))
    hist(purity, xlab = "purity", cex.lab = 1.5, cex.axis = 1.5, col = "black", main = "")
    abline(v = purity_cutoff, col = "red")
  dev.off()  
  
  dz_tissue <-dz_tissue[, purity > purity_cutoff & !is.na(purity)]
  dz_tissue_phenotype <- dz_tissue_phenotype %>% 
    filter(sample %in% colnames(dz_tissue)) #maintain consistency
  rm(sample_purity)
}

####Find and Annotate cancer mutation subtype patient samples####
if (mutation_gene != "" & gdc_project_id != ""){
  print("querying GDC for gene mutation samples")
  gdc_mutated_tumor_pt <-queryGDC(id_mapping_gene_ensembl(mutation_gene), 
                                  gdc_project_id) %>% filter(gene == 1,project==1)
  #Annotating dz_tissue_phenotype with conditions
  
  #TODO consider setting levels here maybe it will set downstream
  dz_tissue_phenotype <- dz_tissue_phenotype %>% mutate(condition = factor(ifelse(patient_id %in% 
                                                           gdc_mutated_tumor_pt$submitter_id, 
                                                         paste0(mutation_gene,"-mutant"),
                                                         paste0("non-",mutation_gene,"-tumor"))))
  
  dz_tissue_phenotype$sample_type = "tumor"
  rm(gdc_mutated_tumor_pt)
  dz_tissue_phenotype <- dz_tissue_phenotype %>% 
    filter(sample %in% colnames(dz_tissue)) #maintain consistency

} else if(mutation_gene == "" & gdc_project_id !=""){
  print("querying GDC for tumor samples")
  gdc_mutated_tumor_pt <-queryGDC(PROJECT = gdc_project_id)
  dz_tissue_phenotype <- dz_tissue_phenotype %>% 
    filter(patient_id %in% gdc_mutated_tumor_pt$submitter_id)
  dz_tissue <- dz_tissue[,colnames(dz_tissue) %in% dz_tissue_phenotype$sample]
  #Annotating dz_tissue_phenotype with conditions
  dz_tissue_phenotype$sample_type = "tumor"
  dz_tissue_phenotype$condition = 'tumor'
  rm(gdc_mutated_tumor_pt)
  dz_tissue_phenotype <- dz_tissue_phenotype %>% 
    filter(sample %in% colnames(dz_tissue)) #maintain consistency
}else {
  dz_tissue_phenotype$sample_type = "tumor"
}

####Grab Adjacent Normal Tissue Cells####
if(ref_selection == 'TCGA'){
  #all_pheno <- read.csv("~/Documents/GitHub/OCTAD v180104/data/raw/treehouse/TcgaTargetGTEX_phenotype.txt",sep="\t")
  ref_samples_phenotype <- all_pheno %>% filter(toupper(primary.disease.or.tissue) == toupper(dz),X_study == 'TCGA', X_sample_type=='Solid Tissue Normal')
  if(length(ref_samples_phenotype$sample)>0){
    ref_tissue <- dz_expr[,colnames(dz_expr) %in% ref_samples_phenotype$sample]
    ref_samples_phenotype$sample_id <- ref_samples_phenotype$sample 
    ref_samples_phenotype$condition = 'normal'
    ref_samples_phenotype$sample_type = 'normal'
    ref_samples_phenotype <- ref_samples_phenotype %>% filter(sample_id %in% colnames(ref_tissue))
  }else{
    print('No normal adjacent tissue found, switching search to GTEX')
    ref_selection = 'GTEX'}
}

####Get GTEX reference via Compute top correlated cell lines or predefined site####
#can turn on to compute linc data set to use 
#turn off to save computation time
#compute_tissue_cell_cor(colnames(dz_tissue))
#compute_tissue_lincs_cell_cor(colnames(dz_tissue))

#find normal sample expression
#note normal tissue also contains cell line so this code may select the cell line type


if (ref_selection == 'GTEX'){
  if (length(site) == 0|site=="") {
    #test to remove primary site
    #all_pheno <- all_pheno %>% filter(X_primary_site != 'Liver')
    print('determining best ref sample tissue site')
    GTEX_refs <- computeRefTissue(expSet = dz_tissue,
                                  method = selectMethod,
                                  site_selection = siteSelect,
                                  varyingGenes = 5000,
                                  cor_cutoff=ref_tissue_cor_cutoff,
                                  output=T,
                                  visualize=T)
    
    #ref_samples_phenotype <- GTEX_phenotype_cor %>% filter(X_primary_site %in% site, cor > ref_tissue_cor_cutoff)
    if(selectMethod == 'random'|siteSelect=='any'){
      ref_samples_phenotype <- all_pheno %>% filter(sample %in% GTEX_refs)
    }else{
      ref_samples_phenotype <- all_pheno %>% filter(sample %in% GTEX_refs$sample)}
    ref_samples_phenotype$sample_id <- ref_samples_phenotype$sample
    #rm(GTEX_refs)
  }else{
    all_normal <- all_pheno %>% filter(X_study == 'GTEX') 
    ref_samples_phenotype <- all_pheno %>% filter(X_study == 'GTEX',primary.disease.or.tissue==site)
    ref_samples_phenotype$sample_id <- ref_samples_phenotype$sample
  }
}


ref_samples_phenotype$condition = 'normal'
ref_samples_phenotype$sample_type = 'normal'

if (length(ref_samples_phenotype$sample_id) < 4) { stop("few reference tissue samples ") }

ref_tissue <- dz_expr[, colnames(dz_expr) %in%  ref_samples_phenotype$sample_id ] 
#rownames(ref_tissue) <- rownames(dz_expr)
ref_samples_phenotype <- ref_samples_phenotype %>% filter(sample_id %in% colnames(ref_tissue))

####Remove Outlier Ref tissue (opt'l)####
if (remove_outlier_ref == T){
  print("removing outlier ref tissue samples")
  outliers <- detectOutlier(ref_tissue,3,outlierPdf = '/ref_tissue_mds.pdf')
  numOutliers = length(outliers)
  ifelse(numOutliers == 0,
         print('no outliers found'),
         print(paste0(numOutliers," outliers removed")))
  ref_tissue = ref_tissue[, !colnames(ref_tissue) %in% outliers]
  rm(numOutliers,outliers)
  ref_samples_phenotype <- ref_samples_phenotype %>% 
    filter(sample %in% colnames(ref_tissue)) #maintain consistency
}

####Creating Counts DataFrame####
#log2(norm_count+1)
dz_tissue <- 2^dz_tissue - 1
ref_tissue <- 2^ref_tissue -1

#should we reformat RSEM count
#write.csv(dz_tissue, "~/Downloads/dz_tissue.csv")
#txi_dz = tximport("~/Downloads/dz_tissue.csv", type = "rsem")
#sampleTable <- data.frame(condition = factor(rep(c("A", "B"), each = 3)))
#rownames(sampleTable) <- colnames(txi$counts)
#dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)

counts <-cbind(dz_tissue, ref_tissue)
#rownames(counts) <-as.character(dz_expr$sample)

if (contrast_mutants == 1) {
  #to make sure that levels are preserved
  #helpful to interpret design matrix if intercept = normal
  counts_phenotype <- dz_tissue_phenotype %>% 
    select(sample_id,condition) %>% 
    rbind( ref_samples_phenotype %>% select(sample_id,condition))

                           
  counts_phenotype$condition <- factor( counts_phenotype$condition,
                                        levels = c('normal', 
                                          paste0("non-",mutation_gene,"-tumor"),
                                          paste0(mutation_gene,"-mutant")))
  
  
  coldata <- data.frame(sample_id = colnames(counts)) %>% left_join(counts_phenotype)
}else{
  counts_phenotype <- rbind(dz_tissue_phenotype %>% dplyr::select(sample_id,sample_type),
                            ref_samples_phenotype %>% dplyr::select(sample_id,sample_type))
  counts_phenotype$condition <- factor(counts_phenotype$sample_type,levels = c('normal','tumor')) #to keep condition column consistent when fed to model.matrix
  coldata <- data.frame(sample_id = colnames(counts)) %>% left_join(counts_phenotype)
}

#hmm.. lots of odd counts?  filtering genes with very large count may miss signficant genes...we may use ruvseqEmpNorm to normalize counts first 
#counts = counts[rowSums(counts) > 0 & rowMax(counts) < 500000, ]
#detect outliers and normalize counts across multiple studies
#need to run the code manually and inspect plots carefully
#need to choose outliers manually
#source("../code/dz/rna_seq_normalization.R") replaced by the function ruvseqEmpNorm

####Remove Lowly Expressed Genes####
if(normalize_samples == T){
highExpGenes <- remLowExpr(counts,coldata)
write.csv(data.frame(highExpGenes),file=paste0(outputFolder,"highExpGenes.csv"))
counts <- counts[highExpGenes,]

empiricalGenes <- compEmpContGenes(counts,coldata,n_topGenes = n_topGenes)
write.csv(data.frame(empiricalGenes),file=paste0(outputFolder,"computedEmpGenes.csv"))
###RUVg and edgeR method###
set <- newSeqExpressionSet(counts = round(counts),
                           phenoData = data.frame(coldata,row.names= coldata$sample_id))

#library(RColorBrewer)

#colors <- brewer.pal(3,"Set2")
#visualize before UQ norm
# plotRLE(set,outline=F,ylim=c(-4,4),col=colors[set$condition])
# plotPCA(set,col=colors[set$condition],cex=0.2)
#upper quantile normalization

#set <- betweenLaneNormalization(set,which = 'upper')
#using betweenLaneNormalization makes error in RUVg for LGG

set1 <- RUVg(set,empiricalGenes,k=k)
rm(set)
#pData(set1) #stores the normalization of set1

#plotRLE(set1,outline=F,ylim=c(-4,4),col=colors[set$condition])
#plotPCA(set1,col=colors[set$condition],cex=0.2)
}else{set <- newSeqExpressionSet(round(counts),
                                 phenoData = data.frame(coldata,row.names= coldata$sample_id))
      set1 <- set}

##edgeR DE implementation##


#removed as most of this code is in remove lowly expressed genes
# load_edgeR_counts <- function(){
#   x <-DGEList(counts = counts, group = coldata$condition )
#   #RUVg(x)
#   #x <- calcNormFactors(x, method = "TMM")
#   #warning output
#   #In .calcFactorWeighted(obs = x[, i], ref = x[, refColumn],  ... : NaNs produced
#   
#   #cpm_x <- cpm(x)
#   #lcpm_x <- cpm(x, log=TRUE)
#   
#   #remove lowly expressed genes
#   ####Possibly where IDH1 may be deleted?####
#   keep.exprs <- rowSums(cpm_x>1) >= min(table(coldata$condition))
#   x <- x[keep.exprs,, keep.lib.sizes=FALSE]
#   x
# }

# 
# if (normalize_samples == T){
#   print('normalizing tissue samples')
#   counts <-ruvseqEmpNorm(counts, coldata, n_topGenes = 20000, k = 1)
# }else{counts <-round(counts)} #round(counts) is ran in ruvseqEmpNorm function

####Create Dz and Ref tissue counts DGEList for edgeR or Limma####


#unable to run this 
# pdf(paste0(outputFolder, "/tissue_normal_mds.pdf"))
#   #col.group <-coldata$condition
#   #levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
#   #col.group <- as.character(col.group)
#   
#   col.sample <- c("purple","orange")[coldata$condition]
#   
#   plotMDS(x, labels = NULL,  pch = 20, col = col.sample)
#   legend("topleft",fill=c("purple","orange"),legend=levels(coldata$condition))
# dev.off()

#interactive plot
#glMDSPlot(x, labels = colnames(x), groups=x$samples, col = col.group, launch=TRUE)


####DE Method Functions, edgeR, limma, deseq####

if (DE_method == "edgeR"){
  print('computing DE via edgeR')
  if(normalize_samples == T){
    if(k==1){
    design <- model.matrix(~condition + W_1, data=pData(set1))
    }else if(k == 2){
      design <- model.matrix(~condition + W_1 + W_2, data = pData(set1))
    }else if (k == 3){
      design <- model.matrix(~condition + W_1 + W_2 + W_3, data = pData(set1))
    }
  }else{design <- model.matrix(~condition,data=pData(set1))}
  
  dgList <- DGEList(counts=counts(set1),group=set1$condition)
  dgList <- calcNormFactors(dgList, method="TMM") #using upperquartile seems to give issue for LGG
  dgList <- estimateGLMCommonDisp(dgList, design)
  dgList <- estimateGLMTagwiseDisp(dgList, design)
  fit <- glmFit(dgList, design)
  #lrt <- glmLRT(fit)

  #see edgeRUsersGuide section on testing for DE genes for contrast
  #fit <- glmFit(dgList,design = designMat)
  if(contrast_mutants==1){
    lrt <- glmLRT(fit,coef=3) #gene mutant tumor vs. ref tissue
    res <- lrt$table
    
    lrt1 <- glmLRT(fit,coef = 2) #non-gene mutant tumor vs. ref tissue
    res1 <- lrt1$table
    
    lrt2 <- glmLRT(fit,contrast=c(0,-1,1,0)) #gene mutant tumor vs. non-gene mutant tumor
    res2 <- lrt2$table 
    
    colnames(res) <- c("log2FoldChange", "logCPM", "LR", "pvalue")
    res$padj <- p.adjust(res$pvalue)
    colnames(res1) <- c("log2FoldChange", "logCPM", "LR", "pvalue")
    res1$padj <- p.adjust(res1$pvalue)
    colnames(res2) <- c("log2FoldChange", "logCPM", "LR", "pvalue")
    res2$padj <- p.adjust(res2$pvalue)
    res$contrast <- 'mutant vs. ref'
    res1$contrast <- 'non-mutant vs. ref'
    res2$contrast <- 'mutant vs. non-mutant'
    res <- rbind(res,res1,res2)
    #write.csv(res,paste0(outputFolder,"all_genes_differences.csv"))
  }else{
    lrt <- glmLRT(fit,2) 
    #second coefficient otherwise it'll default the W_1 term when normalize is on
    #gene mutant tumor vs. ref tissue
    res <- lrt$table
    colnames(res) <- c("log2FoldChange", "logCPM", "LR", "pvalue")
    res$padj <- p.adjust(res$pvalue)
  }
  
  #deciding for significant DE genes
  #deGenes <- decideTestsDGE(lrt,p=0.005)
  #deGenes <- row.names(lrt)[as.logical(deGenes)]
  #plotSmear(lrt,de.tags = deGenes)
  #abline(h=c(-1,1),col=2)
  
}else if (DE_method == "limma"){
  print('computing DE via limma')
  library('Glimma')
  x <- load_edgeR_counts()
  nsamples <- ncol(x)
  col <- brewer.pal(nsamples, "Paired")
  lcpm <- cpm(x, log=TRUE)
  #plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
  #     main="", xlab="")
  #title(main="B. Filtered data", xlab="Log-cpm")
  #abline(v=0, lty=3)
  # for (i in 2:nsamples){
  #   den <- density(lcpm[,i])
  #   lines(den$x, den$y, col=col[i], lwd=2)
  # }
  
  group <-coldata$condition
  col.group <- group
  levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
  col.group <- as.character(col.group)
  #plotMDS(lcpm, labels=coldata$condition, col=col.group)

  # ----design
  design <- model.matrix(~0 + group)
  colnames(design) <- gsub("group", "", colnames(design))
  
  contr.matrix <- makeContrasts(
    TumorvsNon = tumor - normal, 
    levels = colnames(design))
  
  v <- voom(x, design, plot=F)
  vfit <- lmFit(v, design)
  vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
  efit <- eBayes(vfit)
  #plotSA(efit, main="Final model: Mean−variance trend")
  
  tfit <- treat(vfit, lfc=1)
  dt <- decideTests(tfit)
  summary(dt)
  
  tumorvsnormal <- topTreat(tfit, coef=1, n=Inf)
  tumorvsnormal <- tumorvsnormal[order(abs(tumorvsnormal$logFC), decreasing = T),]
  tumorvsnormal.topgenes <- rownames(tumorvsnormal[1:50,])
  'mycol <- colorpanel(1000,"blue","white","red")
  pdf( paste0(outputFolder, "/limma_sig.pdf"))
    heatmap.2(v$E[tumorvsnormal.topgenes,], scale="row",
              labRow=tumorvsnormal.topgenes, labCol=group, 
              col=mycol, trace="none", density.info="none", 
              margin=c(8,6), lhei=c(2,10), dendrogram="column")
  dev.off()'
    
  'df <- as.data.frame(coldata[,c("condition")])
  rownames(df) = (coldata$sample)
  colnames(df) = c("type")
  pheatmap(v$E[tumorvsnormal.topgenes,], cluster_rows=T, show_rownames=T,show_colnames=F, scale="row",col=mycol,
           cluster_cols=T, annotation_col=df, file= paste0(outputFolder, "/limma_sig.pdf"))
  '
  
  res <-tumorvsnormal
  colnames(res) <-c("log2FoldChange", "AveExpr", "t", "pvalue", "padj")
}else{
  library('DESeq2')
  print('computing DE via DESeq')
  dds <- DESeqDataSetFromMatrix(countData = round(counts),
                                colData = coldata,
                                design= ~ condition)
  gc()
  if (parallel_cores > 1){
    dds <- DESeq(dds, parallel = T)
  }else{
    dds <- DESeq(dds)
  }
  
  save(dds,file= paste0(outputFolder, "/dds", ".RData"))
  rnms <- resultsNames(dds)
  
  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:50]
  df <- as.data.frame(colData(dds)[,c("condition")])
  rownames(df) <-rownames(colData(dds))
  colnames(df) <-c("type")
  ntd <- normTransform(dds)
  pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
           cluster_cols=FALSE, annotation_col=df, file= paste0(outputFolder, "/deseq_cluster_by_counts.pdf"))
  
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
           col=colors, file= paste0(outputFolder, "/cluster_by_co_correlation.pdf"))
  '
  
  res <- results(dds, contrast=c("condition","tumor","normal"))
}
  
#plotMA(res, ylim=c(-2,2))
#plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

#write.csv(res, paste0(outputFolder, "/dz_sig_genes_all_",DE_method, ".csv")  )

createSignatureMappings <- function(DE_results){
  mapping <- read.csv(paste0(dataFolder,'raw/gene_info_hs.csv'))
  all_results <- DE_results %>% 
    mutate(Symbol = row.names(DE_results)) %>% 
    left_join(mapping)
  write.csv(all_results,paste0(outputFolder,"all_genes_differences.csv"))
  dz_signature <- DE_results %>% 
    mutate(Symbol = row.names(DE_results)) %>% 
    left_join(mapping)
  dz_signature <-dz_signature[order(dz_signature$log2FoldChange), ]
  dz_signature <-dz_signature[abs(dz_signature$log2FoldChange) > dz_fc_threshold & dz_signature$padj < dz_p_threshold,] 
  #& !is.na(dz_signature$Symbol) & !is.na(dz_signature$padj), ]
  dz_signature$value <-dz_signature$log2FoldChange
  dz_signature$up_down <-ifelse(dz_signature$value > 0, "up", "down")
  dz_signature
}

if (contrast_mutants == 1) {
  dz_signature <- createSignatureMappings(res)
  dz_signature$contrast <- 'mutant vs normal'
  write.csv(dz_signature,paste0(outputFolder,"/dz_sig_genes_mutantTumor_vs_normal_",DE_method,".csv"))
  dz_signature1 <- createSignatureMappings(res1)
  write.csv(dz_signature1, paste0(outputFolder, "/dz_sig_genes_nonmutantTumor_vs_normal_", DE_method,".csv"))
  dz_signature1$contrast <- 'non-mutant tumor vs normal'
  dz_signature2 <- createSignatureMappings(res2)
  write.csv(dz_signature2, paste0(outputFolder, "/dz_sig_genes_mutant_vs_nonmutant", DE_method,".csv")  )
  dz_signature2$contrast <- 'mutant vs nonmutant tumor'
  write.csv(rbind(dz_signature,dz_signature1,dz_signature2), 
            paste0(outputFolder, "/dz_sig_genes_", DE_method,".csv"))
  
}else{
  dz_signature <- createSignatureMappings(res)
  write.csv(dz_signature, paste0(outputFolder, "/dz_sig_genes_", DE_method,".csv"))
}

counts_phenotype <- all_pheno %>% filter(sample %in% counts_phenotype$sample_id)
write.csv(counts_phenotype,paste0(outputFolder,"counts_phenotype.csv"))

#doesn't work
#visualize sig genes using lcpm
# df <- as.data.frame(coldata[,c("condition")])
# rownames(df) <-(coldata$sample)
# colnames(df) <-c("type")
# mycol <- colorpanel(1000,"blue","white","red")
# pheatmap(lcpm[as.character(dz_signature$Symbol),], cluster_rows=T, show_rownames=F,show_colnames=F, scale="row",col=mycol,
#          cluster_cols=T, annotation_col=df, file= paste0(outputFolder, "/sig_genes_lcpm.pdf"))

print('Complete dz signature')
#comparing EdgeR, limma, deseq

#In HCC, limma missed the most well-known biomarker "AFP", deseq is super slow, edgeR relatively better in terms of accuracy and speed
#clear environment
#rm(list=ls())
