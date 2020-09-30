#' @export
diffExp <- function(case_id='',control_id='',source='octad',file='octad.counts.and.tpm.h5',
                    normalize_samples=TRUE,
                    k=1,
					expSet=NULL,
                    n_topGenes=500,
					          n_varGenes=NULL,
                    DE_method='edgeR',
                    parallel_cores = 2,
					output=TRUE,
					outputFolder='', annotate=TRUE){
#  require(dplyr)
#  require(RUVSeq)
#  require(edgeR)
if(!missing(n_varGenes)){
  n_topGenes=n_varGenes
}

if(missing(case_id)|missing(control_id)){
stop('Case ids and/or control ids vector input not found')
}


if(source=='octad.whole'){

print(paste('loading whole octad expression data for',length(c(case_id,control_id)),'samples',sep=' '))

#case.load.start <- Sys.time()
#setwd(base.folder)
#rhdf5_file <- paste0(, "octad.h5")
#transcripts=row.names(octad.db::octad_counts)
#saamples=colnames(octad.db::octad_counts)

#expSet=octad.db::octad_counts[row.names()]
  



transcripts = as.character(rhdf5::h5read(file, "meta/transcripts"))
samples = as.character(rhdf5::h5read(file, "meta/samples"))
case_counts = rhdf5::h5read(file, "data/count", index=list(1:length(transcripts), which(samples %in% case_id)))
colnames(case_counts) = samples[samples %in% case_id]
rownames(case_counts) = transcripts
case_id = samples[samples %in% case_id]
normal_counts = rhdf5::h5read(file, "data/count", index=list(1:length(transcripts), which(samples %in% control_id)))
colnames(normal_counts) = samples[samples %in% control_id]
rownames(normal_counts) = transcripts
control_id = samples[samples %in% control_id]
H5close()

expSet = cbind(normal_counts, case_counts)
rm(normal_counts) # free up some memory
rm(case_counts)

}else if(source=='octad.small'){
	print('loading small octad set containing only expression for 978 LINCS genes')
#case.load.start <- Sys.time()
#setwd(base.folder)
#rhdf5_file <- paste0(, "octad.h5")
#transcripts = as.character(rhdf5::h5read(file, "meta/transcripts"))
#samples = as.character(rhdf5::h5read(file, "meta/samples"))
#case_counts = rhdf5::h5read(file, "data/count", index=list(1:length(transcripts), which(samples %in% case_id)))
#colnames(case_counts) = samples[samples %in% case_id]
#rownames(case_counts) = transcripts
#case_id = samples[samples %in% case_id]
#normal_counts = rhdf5::h5read(file, "data/count", index=list(1:length(transcripts), which(samples %in% control_id)))
#colnames(normal_counts) = samples[samples %in% control_id]
#rownames(normal_counts) = transcripts
#control_id = samples[samples %in% control_id]
#H5close()
expSet=octad.db::octad.LINCS.counts[,c(case_id,control_id)]
#expSet = cbind(normal_counts, case_counts)
#rm(normal_counts) # free up some memory
#rm(case_counts)

}else if(source!='octad.small'&source!='octad.small'&missing(expSet)){
stop('Expression data not sourced, please, modify expSet option')

}




  counts_phenotype <- rbind(data.frame(sample = case_id,sample_type = 'case'),
                            data.frame(sample = control_id, sample_type = 'control'))
  counts_phenotype = counts_phenotype[counts_phenotype$sample %in% colnames(expSet),]
  
  counts = expSet[,as.character(counts_phenotype$sample)]
  counts = 2^counts - 1 #unlog the counts it was log(2x + 1) in dz.expr.log2.readCounts
  counts_phenotype$sample = as.character(counts_phenotype$sample)
  counts_phenotype$sample_type = factor(counts_phenotype$sample_type, levels = c("control", "case"))
  
  #remove lowly expressed transcripts
  highExpGenes <- remLowExpr(counts,counts_phenotype)
  counts = counts[highExpGenes,]
  
  set <- EDASeq::newSeqExpressionSet(round(counts),
                             phenoData = data.frame(counts_phenotype,row.names=counts_phenotype$sample))
  
  #normalize samples using RUVSeq
  if (normalize_samples == T){
    #compute empirical genes
    
    design <- stats::model.matrix(~ sample_type, data = pData(set))
    y <- edgeR::DGEList(counts=counts(set), group =  counts_phenotype$sample)
    y <- edgeR::calcNormFactors(y, method="TMM") #upperquartile generate Inf in the LGG case
    y <- edgeR::estimateGLMCommonDisp(y, design)
    y <- edgeR::estimateGLMTagwiseDisp(y, design)
    fit <- edgeR::glmFit(y, design)
    lrt <- edgeR::glmLRT(fit,2) #defaults to compare case control
    
    top <- edgeR::topTags(lrt, n=nrow(set))$table
    i = which(!(rownames(set) %in% rownames(top)[1:min(n_topGenes,dim(top)[1])]))
    empirical <- rownames(set)[i]
    stopifnot(length(empirical)>0)
    write.csv(data.frame(empirical),file=paste0(outputFolder,"computedEmpGenes.csv"))
    set1 <- RUVSeq::RUVg(set,empirical,k=k)
  }
  
  if(DE_method=='DESeq2'){
#    library('DESeq2')
    print('computing DE via DESeq')
    row.names(counts_phenotype) = counts_phenotype$sample
    coldata = counts_phenotype
    
    if (normalize_samples == T){
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts(set1),
                                    colData = pData(set1),
                                    design= ~ sample_type + W_1)
    }else{
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(counts),
                                    colData = coldata,
                                    design= ~ sample_type)
    }
    gc()
    
    if (parallel_cores > 1){
      dds <- DESeq2::DESeq(dds, parallel = T)
    }else{
      dds <- DESeq2::DESeq(dds)
    }
    
    #save(dds,file= paste0(outputFolder, "/dds", ".RData"))
    rnms <- DESeq2::resultsNames(dds)
    
    resRaw <- DESeq2::results(dds, contrast=c("sample_type","case","control"))
    res = data.frame(resRaw)
    res$identifier <- row.names(res)
    res = res %>% select(identifier,everything())
#    return(res)
    
  }else if(DE_method=='edgeR'){
    print('computing DE via edgeR')
    
    #construct model matrix based on whether there was normalization ran
    if(normalize_samples == T){
      if(k==1){
        design <- stats::model.matrix(~sample_type + W_1, data=pData(set1))
      }else if(k == 2){
        design <- stats::model.matrix(~sample_type + W_1 + W_2, data = pData(set1))
      }else if (k == 3){
        design <- stats::model.matrix(~sample_type + W_1 + W_2 + W_3, data = pData(set1))
      }
      dgList <- edgeR::DGEList(counts=counts(set1),group=set1$sample_type)
      
    }else{
      design <- stats::model.matrix(~ sample_type,data=pData(set))
      dgList <- edgeR::DGEList(counts=counts(set),group=set$sample_type)
    }
    dgList <- edgeR::calcNormFactors(dgList, method="TMM") #using upperquartile seems to give issue for LGG
    dgList <- edgeR::estimateGLMCommonDisp(dgList, design)
    dgList <- edgeR::estimateGLMTagwiseDisp(dgList, design)
    fit <- edgeR::glmFit(dgList, design)
    #see edgeRUsersGuide section on testing for DE genes for contrast
    lrt <- edgeR::glmLRT(fit,2) 
    #second coefficient otherwise it'll default the W_1 term when normalize is on
    res <- lrt$table
    colnames(res) <- c("log2FoldChange", "logCPM", "LR", "pvalue")
    res$padj <- p.adjust(res$pvalue)
    res$identifier <- row.names(res)
    res = res %>% select(identifier,everything())
#    return(res)
  }else if(DE_method=='limma'){
    #according to https://support.bioconductor.org/p/86461/, LIMMA + VOOM will not use normalized data
    print('computing DE via limma')
#    require('Glimma')
    x <- counts
    nsamples <- ncol(x)
    lcpm <- edgeR::cpm(x, log=TRUE)
    
    group <-counts_phenotype$sample_type
    
    ## ----design-----------------------------------------------------------------------------
    design <- model.matrix(~0 + group)
    colnames(design) <- gsub("group", "", colnames(design))
    
    contr.matrix <- limma::makeContrasts(
      TumorvsNon = case - control, 
      levels = colnames(design))
    
    v <- limma::voom(x, design, plot=F)
    vfit <- limma::lmFit(v, design)
    vfit <- limma::contrasts.fit(vfit, contrasts=contr.matrix)
    efit <- limma::eBayes(vfit)
    #plotSA(efit, main="Final model: Mean−variance trend")
    
    tfit <- limma::treat(vfit, lfc=1) #not_sure
    dt <- limma::decideTests(tfit)
    summary(dt)
    
    tumorvsnormal <- limma::topTreat(tfit, coef=1, n=Inf)
    tumorvsnormal <- tumorvsnormal[order(abs(tumorvsnormal$logFC), decreasing = T),]
    #tumorvsnormal.topgenes <- rownames(tumorvsnormal[1:50,])
    
    
    res <-tumorvsnormal
    colnames(res) <-c("log2FoldChange", "AveExpr", "t", "pvalue", "padj")
    res$identifier = row.names(res)
#    return(res)
  }
#load(merged_gene_info)
if(annotate==TRUE){
	merged_gene_info=octad.db::merged_gene_info
    merged_gene_info$ensembl<-as.vector(merged_gene_info$ensembl)
    merged_gene_info$V1=NULL #modify it when will have time
    merged_gene_info$Symbol=merged_gene_info$gene #modify it when will have time
    merged_gene_info$gene=NULL
    res <- left_join(res, merged_gene_info, by=c('identifier'='ensembl'))
}
return(res)
}