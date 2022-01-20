#' @export
#' @importFrom Rfast colRanks 
#' @importFrom dplyr summarise id desc
#' @importFrom plotly group_by
#' @import octad.db
#' @importFrom ExperimentHub ExperimentHub


#### runsRGES #######
runsRGES <- function(dz_signature=NULL,choose_fda_drugs = FALSE,max_gene_size=500, 
                                         cells=NULL,outputFolder=NULL,weight_cell_line=NULL,permutations=10000){

if(missing(dz_signature)){
	stop('Disease signature input not found')
}
if(is.null(dz_signature$Symbol)|is.null(dz_signature$log2FoldChange)){
	stop('Either Symbol or log2FoldChange collumn in Disease signature is missing')
}

eh=ExperimentHub()
lincs_sig_info=eh[["EH7270"]] #bioconductor addition
getsRGES <- function(RGES, cor, pert_dose, pert_time, diff, max_cor){        
	sRGES <- RGES
	pert_time <- ifelse(pert_time < 24, "short", "long")
	pert_dose <- ifelse(pert_dose < 10, "low", "high")
	if (pert_time == "short" & pert_dose == "low"){
		sRGES <- sRGES + diff[4]
        }
	if (pert_dose ==    "low" & pert_time == "long"){
		sRGES <- sRGES + diff[2]
        }
	if (pert_dose ==    "high" & pert_time == "short"){
		sRGES <- sRGES + diff[1]
        }
	return(sRGES * cor/max_cor) 
    }
    
cmap_score_ultimate=function (sig_up, sig_down, drug_signature){
	num_genes <- length(drug_signature)
	ks_up <- 0
	ks_down <- 0
	connectivity_score <- 0
    drug_signature=rank(drug_signature)
    up_tags_rank=drug_signature[as.vector(sig_up)]
    down_tags_rank=drug_signature[as.vector(sig_down)]
	up_tags_position=sort(up_tags_rank)
	down_tags_position=sort(down_tags_rank)
	num_tags_up <- length(up_tags_position)
	num_tags_down <- length(down_tags_position)
	if (num_tags_up > 1) {
		a_up <- 0
		b_up <- 0
		a_up <- max(sapply(seq_len(num_tags_up), function(j) {
			j/num_tags_up - up_tags_position[j]/num_genes
			}))
		b_up <- max(sapply(seq_len(num_tags_up), function(j) {
			up_tags_position[j]/num_genes - (j - 1)/num_tags_up
			}))
	if (a_up > b_up) {
		ks_up <- a_up
		}else{
			ks_up <- -b_up
	}}
	else {
		ks_up <- 0
        }
if (num_tags_down > 1) {
		a_down <- 0
		b_down <- 0
		a_down <- max(sapply(seq_len(num_tags_down), function(j) {
			j/num_tags_down - down_tags_position[j]/num_genes
		}))
		b_down <- max(sapply(seq_len(num_tags_down), function(j) {
			down_tags_position[j]/num_genes - (j - 1)/num_tags_down
		}))
if (a_down > b_down) {
		ks_down <- a_down
	}
	else {
		ks_down <- -b_down
	}}
	else {
		ks_down <- 0
	}
if (ks_up == 0 & ks_down != 0) {
	connectivity_score <- -ks_down
	}
	else if (ks_up != 0 & ks_down == 0) {
		connectivity_score <- ks_up
        }
	else if (sum(sign(c(ks_down, ks_up))) == 0) {
		connectivity_score <- ks_up - ks_down
        }
	else {
		connectivity_score <- ks_up - ks_down
        }
        return(connectivity_score)
}
  
if (!missing(outputFolder)) {
	if(!dir.exists(outputFolder)){
        dir.create(outputFolder)
        }
    }

if(!missing(cells)){		
#lincs_sig_info=ExperimentHub()[["EH7270"]]
	lincs_sig_info$cell_id = toupper(lincs_sig_info$cell_id)
	#lincs_sig_info$cell_id = toupper(octad.db::lincs_sig_info$cell_id) #bioconductor replace		
	lincs_sig_info = subset(lincs_sig_info,cell_id %in% cells)
	#lincs_sig_info = subset(octad.db::lincs_sig_info,cell_id %in% cells) #bioconductor replace		
    }else if(!missing(cells)&nrow(lincs_sig_info)==0){
        stop('Wrong cell line name. List of possible cell lines is available via command unique(octad.db::lincs_sig_info$cell_id)')
    }else if(missing(cells)){
        cells='' #plug for older code version
    }		

lincs_signatures=eh[["EH7271"]]     

if (choose_fda_drugs) {	
	fda_drugs=eh[["EH7269"]] 
    #fda_drugs=octad.db::fda_drugs #bioconductor replace	
	lincs_sig_info_FDA <- subset(lincs_sig_info, id %in% colnames(lincs_signatures) & tolower(pert_iname) %in% tolower(fda_drugs$pert_iname))
	FDAdf <- select(lincs_sig_info_FDA, pert_id, pert_iname)
	FDAdf <- unique(FDAdf[,seq_len(2)])
	write.csv(FDAdf,file = paste0(outputFolder,"FDA_approved_drugs.csv"),row.names = FALSE)
	lincs_sig_info <- subset(lincs_sig_info,id %in% colnames(lincs_signatures))
    }else{
        lincs_sig_info <- subset(lincs_sig_info, id %in% colnames(lincs_signatures))
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
dz_signature <- subset(dz_signature,Symbol %in% gene.list)
dz_genes_up <- subset(dz_signature,log2FoldChange>0)
dz_genes_up=dz_genes_up[order(dz_genes_up$log2FoldChange,decreasing=TRUE),]
dz_genes_down <- subset(dz_signature,log2FoldChange<0)
dz_genes_down=dz_genes_down[order(dz_genes_down$log2FoldChange,decreasing=TRUE),]
###############
    
    
#compute RGES
#caps gene selection to max gene size 
if (nrow(dz_genes_up) > max_gene_size){
	dz_genes_up <- dz_genes_up %>% head(max_gene_size)
    }
if (nrow(dz_genes_down) > max_gene_size){
	dz_genes_down <- dz_genes_down %>% head(max_gene_size)
    }
    
write.csv(rbind(dz_genes_up, dz_genes_down),    dz_sig_output_path)
      
#parallel=FALSE
#if(!parallel){ in case of parallel implementation for later release
dz_cmap_scores <- NULL
	#count <- 0
cmap_exp_sig <- Rfast::colRanks(-1 * lincs_signatures, method = "max")    
names.list <- list(rownames(lincs_signatures),colnames(lincs_signatures))
dimnames(cmap_exp_sig) <- names.list
cat(paste('Started sRGES computation. Average computation time ~1-3mins.'),'\n')
start_time=Sys.time()   
pb <- txtProgressBar(min = 1, max = permutations, style = 3) #set progressbar
i=0
time_vector=0
   
cmap_exp_signature=data.frame(ids = gene.list, rank = cmap_exp_sig[,as.vector(sig.ids)])
cmap_exp_sig=as.data.frame(cmap_exp_sig)
cmap_exp_sig$ids=NULL
dz_cmap_scores=apply(cmap_exp_sig[as.vector(sig.ids)],2,
            FUN=function(x) cmap_score_ultimate(dz_genes_up$Symbol,dz_genes_down$Symbol,drug_signature=x))


#random scores
lincs_signatures=eh[["EH7271"]] 
random_sig_ids <- sample(colnames(lincs_signatures),permutations,replace=TRUE)
random_cmap_scores <- NULL

cmap_exp_signature=as.data.frame(Rfast::colRanks(-1 * lincs_signatures[, as.character(random_sig_ids)], method = "max"))
rm(lincs_signatures) #free some memory

random_cmap_scores=apply(cmap_exp_signature,2,
	FUN=function(x) cmap_score_ultimate(
	sample(1:length(dz_genes_up$Symbol),replace=TRUE),
	sample(1:length(dz_genes_down$Symbol),replace=TRUE),
	drug_signature=x))

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
        
lincs_drug_prediction_subset <- subset(lincs_drug_prediction,    pert_dose > 0 & pert_time %in% c(6, 24))
#pairs that share the same drug and cell id
lincs_drug_prediction_pairs <- merge(lincs_drug_prediction_subset, lincs_drug_prediction_subset, by=c("pert_iname", "cell_id")) 
#x is the reference
lincs_drug_prediction_pairs <- subset(lincs_drug_prediction_pairs, id.x != id.y & pert_time.x == 24 & pert_dose.x == 10) #, select <- c("cmap_score.x", "cmap_score.y", "pert_dose.y", "pert_time.y"))
        
#difference of RGES to the reference 
lincs_drug_prediction_pairs$cmap_diff <- lincs_drug_prediction_pairs$cmap_score.x - lincs_drug_prediction_pairs$cmap_score.y
lincs_drug_prediction_pairs$dose <- round(log(lincs_drug_prediction_pairs$pert_dose.y, 2), 1)
#cat(paste('mem check 13'),'\n')       
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

if (choose_fda_drugs){
	pred_merged_drug <- merge(fda_drugs, pred_merged, by = "pert_iname")
	pred_merged_drug <- pred_merged_drug[order(pred_merged_drug$sRGES), ]
	write.csv(pred_merged_drug, sRGES_output_path_drug)
	}
cat('\n',paste('Finished computations in', round(Sys.time()-start_time,2),units(Sys.time()-start_time),',writing output'),'\n')
return(pred_merged)
#    } end of do.parallel
}
