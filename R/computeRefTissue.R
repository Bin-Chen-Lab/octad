#' @export
#' @importFrom magrittr %>%	
computeRefTissue <- function(case_id = NULL,
								adjacent=FALSE,
                             source='octad',n_varGenes = 500,
                             method='varGenes', #random 
							 expSet=NULL,
                             #site_selection = 'top', 
                             #top site or any cor cutoff>quantile,
                             #all will select samples from 90th percentile
                             control_size = length(case_id),
							 outputFolder='',
                             cor_cutoff='0%', #greater or equal than the cutoff 
                             output=TRUE){
	 
#  require(dplyr)
require(octad.db)
if(missing(case_id)){
stop('Case ids vector input not found')
}
if(source=='octad'){
expSet=octad.db::EncoderDF
case_id=case_id[case_id %in% colnames(expSet)] 
#if we pick adjacent, filter them out
if(adjacent==TRUE){

adjacent_ids = as.vector((octad.db::phenoDF %>% dplyr::filter(sample.type == "adjacent"))$sample.id)
normal_id = as.vector((octad.db::phenoDF %>% dplyr::filter(sample.type == "normal"))$sample.id)
normal_id = c(adjacent_ids,normal_id)

}else{
 #load autoencoder dataset as ExpSet
#load('data/EncoderDF.rda')
normal_id = as.vector((octad.db::phenoDF %>% dplyr::filter(sample.type == "normal"))$sample.id)

}


}else if(source!='octad'&missing(expSet)){
stop('expSet is not supported')
}else if(source!='octad'){  
normal_id=colnames(expSet)[!colnames(expSet)%in%case_id]
}

normal_id=normal_id[normal_id%in%colnames(expSet)]


#rm(phenoDF)

#rm('data/EncoderDF.rda')

  if(method == 'random'){
    GTEXid <- sample(normal_id,size = control_size)
    return(GTEXid)
  }else if(method == 'varGenes'){
    expSet_normal <- expSet[,as.vector(normal_id)]
    expSet_case <- expSet[,as.vector(case_id)]
    #varGenes look at the top varying genes (IQR) within normal tissue expression and varies them to the case tissues
    iqr_gene <-apply(expSet_normal, 1, stats::IQR) #get the IQR per gene
    varying_genes <-order(iqr_gene, decreasing=TRUE)[1:min(n_varGenes,length(iqr_gene))]
    
    #get the correlation matrix for each normal id and each case id
    normal_dz_cor <-cor(expSet_normal[varying_genes, ], expSet_case[varying_genes, ], method = "spearman")
    normal_dz_cor_each <-apply(normal_dz_cor, 1, median) #getting the median correlation btw each normal tissue to the case overall
    normal_dz_cor_eachDF = data.frame(cor=sort(normal_dz_cor_each, decreasing=TRUE)) %>% 
      dplyr::mutate(sample.id = row.names(.)) %>% dplyr::select(sample.id,cor)
    cutoff = stats::quantile(normal_dz_cor_eachDF$cor,probs=seq(0,1,0.05),na.rm=TRUE)[cor_cutoff]
    GTEXid <- (normal_dz_cor_eachDF %>% 
                 dplyr::arrange(desc(cor)) %>% 
                 dplyr::filter(cor>=cutoff))$sample.id 
    GTEXid <- GTEXid[1:min(control_size,length(GTEXid))]
    
    if(output==TRUE){
		if(nchar(outputFolder)>0){
		  if (!dir.exists(outputFolder)) {
			dir.create(outputFolder)
			}
		tryCatch(write.csv(normal_dz_cor,file = paste0(outputFolder,'/case_normal_corMatrix.csv')),
               error = function(c) "failed to write case normal cor matrix csv. Try checking if your outputFolder string is correct or exists")
      
		tryCatch(write.csv(normal_dz_cor_eachDF,row.names = F, paste0(outputFolder, "/case_normal_median_cor.csv")),
               error = function(c) "failed to write case normal median correlation csv. Try checking if your outputFolder string is correct or exists")
}else{
    tryCatch(write.csv(normal_dz_cor,file = paste0(outputFolder,'case_normal_corMatrix.csv')),
               error = function(c) "failed to write case normal cor matrix csv. Try checking if your outputFolder string is correct or exists")
      
    tryCatch(write.csv(normal_dz_cor_eachDF,row.names = F, paste0(outputFolder, "case_normal_median_cor.csv")),
               error = function(c) "failed to write case normal median correlation csv. Try checking if your outputFolder string is correct or exists")
    }}
    return(GTEXid)
  }
    if(output==TRUE){
		if(nchar(outputFolder)>0){
		  if (!dir.exists(outputFolder)) {
			dir.create(outputFolder)
			}
		tryCatch(write.csv(normal_dz_cor,file = paste0(outputFolder,'/case_normal_corMatrix.csv')),
               error = function(c) "failed to write case normal cor matrix csv. Try checking if your outputFolder string is correct or exists")
      
		tryCatch(write.csv(normal_dz_cor_eachDF,row.names = F, paste0(outputFolder, "/case_normal_median_cor.csv")),
               error = function(c) "failed to write case normal median correlation csv. Try checking if your outputFolder string is correct or exists")
}else{
    tryCatch(write.csv(normal_dz_cor,file = paste0(outputFolder,'case_normal_corMatrix.csv')),
               error = function(c) "failed to write case normal cor matrix csv. Try checking if your outputFolder string is correct or exists")
      
    tryCatch(write.csv(normal_dz_cor_eachDF,row.names = F, paste0(outputFolder, "case_normal_median_cor.csv")),
               error = function(c) "failed to write case normal median correlation csv. Try checking if your outputFolder string is correct or exists")
    }}
    return(GTEXid)
}
