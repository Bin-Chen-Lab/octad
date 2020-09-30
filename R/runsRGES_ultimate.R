#' @export
#### runsRGES #######
runsRGES <- function(dz_signature=NULL,choose_fda_drugs = F,max_gene_size=500, 
                     cells=NULL,outputFolder=NULL,weight_cell_line=NULL,permutations=10000
){
require(octad.db)
  if(missing(dz_signature)){
stop('Disease signature input not found')
}

if(is.null(dz_signature$Symbol)|is.null(dz_signature$log2FoldChange)){
stop('Either Symbol or log2FoldChange collumn in Disease signature is missing')
}

  if (!missing(outputFolder)) {
	if(!dir.exists(outputFolder)){
		dir.create(outputFolder)
		}
  }
  
  #  require("dplyr")
  #  require("ggplot2")
  #  require(data.table)
  #RGES
  #Compute RGES by calculating reversal gene expression from LINCS1000
  ####DataFrames####
  #lincs_signature : dataframe for lincs pertubation data 
  #row id : landmark genes 978 of them
  #need metadata to change them to Symbol
  #col id : experiment id
  #need metadata to parse them
  #lincs_sig_info : experiment info for lincs pertubation data
  #row id corresponds to lincs_signature column id
  #fda_drugs : fda drug data
  
  ####required input####
  #outputFolder : folder to output the drug resutls
  #dataFolder : root folder with lincs input data
  #dz_signature : disease signature genes dataframe
  #required columns Symbol : this must be an UPPERCASE gene symbol needed to join with RGES landmark genes
  #must change this if they require another identifier
  
  ####parameters####
  #store these parameters somewhere if you are doing multiple runs for comparison
  
  # default parameters
  # landmark = 1
  # choose_fda_drugs = F
  # max_gene_size = 100
  # weight_cell_line = F
  #load dz_signature
  #load output folder
  #load LINCS drug gene expression profiles
  # load("data/lincs_sig_info")
  #  lincs_sig_info <- data.table::fread("data/lincs_sig_info.csv", stringsAsFactors = TRUE)

  if(!missing(cells)){
    lincs_sig_info$cell_id = toupper(octad.db::lincs_sig_info$cell_id)
    lincs_sig_info = octad.db::lincs_sig_info %>% filter(cell_id %in% cells)
  }else if(!missing(cells)&nrow(lincs_sig_info)==0){
    stop('Wrong cell line name. List of possible cell lines is available via command unique(octad.db::lincs_sig_info$cell_id)')
  }else if(missing(cells)){
    cells='' #plug for older code version
  }
#  landmark = 1
 # if (landmark == 1){
    #there are a few lincs dataframes in the datafolder but for some reason only this one has gene symbols...
    #load(paste0('data/lincs_signatures.rda'))
    lincs_signatures=octad.db::lincs_signatures
#  }else{
    #don't bother
 #   load(paste0(data,"lincs_signatures_cmpd_landmark_GSE92742.RData"))
 # }
  
  if (choose_fda_drugs) {
    #fda_drugs = read.csv("data/repurposing_drugs_20170327.txt", stringsAsFactors = F,sep='\t')
    #load('data/fda_drugs.rda')
    fda_drugs=octad.db::fda_drugs
    lincs_sig_info_FDA <- subset(lincs_sig_info, id %in% colnames(lincs_signatures) & tolower(pert_iname) %in% tolower(fda_drugs$pert_iname))
    FDAdf <- select(lincs_sig_info_FDA, pert_id, pert_iname)
    FDAdf <- unique(FDAdf[,1:2])
    write.csv(FDAdf,file = paste0(outputFolder,"FDA_approved_drugs.csv"),row.names = F)
    lincs_sig_info <- lincs_sig_info %>% filter(id %in% colnames(lincs_signatures))
  }else{
    lincs_sig_info <- lincs_sig_info %>% filter(id %in% colnames(lincs_signatures))
  }
  
  
  #write paths
  output_path <- paste0("all_",paste(cells,collapse='_'),"_lincs_score.csv")
  sRGES_output_path <- paste0("sRGES",paste(cells,collapse='_'),".csv")
  sRGES_output_path_drug <- paste0("sRGES_FDAapproveddrugs.csv")
  dz_sig_output_path <- paste0("dz_sig_used.csv")
  
  #remove duplicate instances
  lincs_sig_info <- lincs_sig_info[!duplicated(lincs_sig_info$id),]
  sig.ids <- lincs_sig_info$id
  
  ####compute RGES####
  gene.list <- toupper(rownames(lincs_signatures))
  
  dz_signature <- dz_signature %>% filter(Symbol %in% gene.list)
  dz_genes_up <- dz_signature %>% filter(log2FoldChange>0) %>% arrange(desc(log2FoldChange))
  dz_genes_down <- dz_signature %>% filter(log2FoldChange<0) %>% arrange(log2FoldChange)
  ###############
  
  
  #compute RGES
  #caps gene selection to max gene size 
  if (nrow(dz_genes_up) > max_gene_size){
    dz_genes_up <- dz_genes_up %>% head(max_gene_size)
  }
  if (nrow(dz_genes_down) > max_gene_size){
    #dz_genes_down <- data.frame(GeneID=dz_genes_down[1:max_gene_size,]) %>% left_join(dz_signature,by = 'GeneID')
    dz_genes_down <- dz_genes_down %>% head(max_gene_size)
  }
  
  write.csv(rbind(dz_genes_up, dz_genes_down),  dz_sig_output_path)
  
  #print(paste('finished loading in', round(Sys.time()-start,2),units(Sys.time()-start)))
  #start=Sys.time()
  
  parallel=F
  if(!parallel){
    #    require(lme4)
    #   require(Rfast)
    dz_cmap_scores <- NULL
    #count <- 0
    
    cmap_exp_sig <- Rfast::colRanks(-1 * lincs_signatures, method = "max")  
    names.list <- list(rownames(lincs_signatures),colnames(lincs_signatures))
    dimnames(cmap_exp_sig) <- names.list
    cat(paste('Started sRGES computation. Average computation time ~1-3mins.'),'\n')
    start_time=Sys.time()  
    ####slow loop#### NO SLOW LOOP ANYMORE! :)
    pb <- txtProgressBar(min = 1, max = permutations, style = 3) #set progressbar
    i=0
    time_vector=0
    
    
#loop_time=Sys.time()
cmap_exp_signature=data.frame(ids = gene.list, rank = cmap_exp_sig[,as.vector(sig.ids)])
cmap_exp_sig=as.data.frame(cmap_exp_sig)
cmap_exp_sig$ids=NULL
dz_cmap_scores=apply(cmap_exp_sig[as.vector(sig.ids)],2,
      FUN=function(x) cmap_score_ultimate(dz_genes_up$Symbol,dz_genes_down$Symbol,drug_signature=x))

#for (exp_id in sig.ids) {
 # i=i+1
 # setTxtProgressBar(pb, i) 
  #count <- count + 1
  #print(count)
  
  #cmap_exp_signature is a dataframe with columns ids: whatever id we're using e.g. Symbol
  #rank which is the reverse rank of the expression of the expression
  # cmap_exp_signature <- data.frame(gene.list,
  #                                 Rfast::colRanks(-1 * lincs_signatures[, as.character(exp_id)],
  #                                          method = "min"))    
  #  colnames(cmap_exp_signature) <- c("ids","rank")
  #runs a function cmap_score_new from drugs_core_functions.R
  #	  setTxtProgressBar(pb, exp_id) 
#  x=Sys.time()
  #  cmap_exp_signature <- data.frame(ids = gene.list, rank = cmap_exp_sig[,exp_id])
#  dz_cmap_scores <- c(dz_cmap_scores, 
#                      cmap_score_new(dz_genes_up$Symbol,dz_genes_down$Symbol,
#                                     cmap_exp_signature[c('ids',paste('rank',exp_id,sep='.'))]))
#  time_vector=c(time_vector,Sys.time()-x)
#}
#Sys.time()-loop_time   
    #print(paste('finished slow loop in', round(Sys.time()-start,2),units(Sys.time()-start)))
    #start=Sys.time()    
    
	
	
    #random scores
#    N_PERMUTATIONS <- 10000 #default 100000

random_sig_ids <- sample(colnames(lincs_signatures),permutations,replace=TRUE)
#count <- 0
random_cmap_scores <- NULL

cmap_exp_signature=as.data.frame(Rfast::colRanks(-1 * lincs_signatures[, as.character(random_sig_ids)], method = "max"))

random_cmap_scores=apply(cmap_exp_signature,2,
                         FUN=function(x) cmap_score_ultimate(
                           sample(1:length(dz_genes_up$Symbol),replace=TRUE),
                           sample(1:length(dz_genes_down$Symbol),replace=TRUE),
                           drug_signature=x))

#random_input_signature_genes <- sample(gene.list, (nrow(dz_genes_up)+nrow(dz_genes_down)))

#rand_dz_gene_up <- data.frame(GeneID=random_input_signature_genes[1:nrow(dz_genes_up)])
#rand_dz_gene_down <- data.frame(GeneID=random_input_signature_genes[(nrow(dz_genes_up)+1):
#length(random_input_signature_genes)])

#random_cmap_scores <- c(random_cmap_scores, cmap_score_new(rand_dz_gene_up,rand_dz_gene_down,cmap_exp_signature))
   #print(paste('ready to start permutations', round(Sys.time()-start,2),units(Sys.time()-start)))
    #start=Sys.time()
#    for (expr_id in random_sig_ids){
 #     i=i+1
 #     setTxtProgressBar(pb, i) 
 #     count <- count + 1
      #print(count)
  #    cmap_exp_signature <- data.frame(gene.list,  
 #                                      rank(-1 * lincs_signatures[, as.character(expr_id)], ties.method="random"))    
  #    colnames(cmap_exp_signature) <- c("ids","rank")
      # cmap_exp_sig <- Rfast::colRanks(-1 * lincs_signatures, method = "min")    
      # names.list <- list(rownames(lincs_sig),colnames(lincs_sig))
      # dimnames(cmap_exp_signature) <- names.list
      
#      random_input_signature_genes <- sample(gene.list, (nrow(dz_genes_up)+nrow(dz_genes_down)))
#      rand_dz_gene_up <- data.frame(GeneID=random_input_signature_genes[1:nrow(dz_genes_up)])
#      rand_dz_gene_down <- data.frame(GeneID=random_input_signature_genes[(nrow(dz_genes_up)+1):length(random_input_signature_genes)])
#      random_cmap_scores <- c(random_cmap_scores, cmap_score_new(rand_dz_gene_up,rand_dz_gene_down,cmap_exp_signature))
#    }
    #print(paste('finished lpermutations', round(Sys.time()-start,2),units(Sys.time()-start)))
    #start=Sys.time()   
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
    if (!missing(weight_cell_line)){
      lincs_cell_line_weight <- weight_cell_line
      pred = merge(lincs_drug_prediction, lincs_cell_line_weight, by.x="cell_id", by.y = 0)
    }else{
      pred <- lincs_drug_prediction
      pred$medcor <- 1
    }
    pred$RGES <- sapply(1:nrow(pred), function(id){getsRGES(pred$cmap_score[id], pred$medcor[id], pred$pert_dose[id], pred$pert_time[id], diff, max(pred$medcor))})
    
    cmpd_freq <- table(pred$pert_iname)
    pred <- subset(pred, pert_iname %in% names(cmpd_freq[cmpd_freq > 0]))
    
    pred_merged <- pred %>% 
      group_by(pert_iname) %>% 
      dplyr::summarise(
        mean = mean(RGES),
        n = length(RGES),
        median = median(RGES),
        sd = sd(RGES))
    pred_merged$sRGES <- pred_merged$mean
    pred_merged <- pred_merged[order(pred_merged$sRGES), ]
    write.csv(pred_merged, sRGES_output_path)
#    return(pred_merged)
    #limit to FDA approved drugs
    if (choose_fda_drugs){
      pred_merged_drug <- merge(fda_drugs, pred_merged, by = "pert_iname")
      pred_merged_drug <- pred_merged_drug[order(pred_merged_drug$sRGES), ]
      write.csv(pred_merged_drug, sRGES_output_path_drug)
    }
    cat('\n',paste('Finished computations in', round(Sys.time()-start_time,2),units(Sys.time()-start_time),',writing output'),'\n')
    #start=Sys.time()
	return(pred_merged)
  }
}
