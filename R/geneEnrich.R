#' @export
geneEnrich=function(gene_list=NULL,database_list=NULL,output=FALSE){
  if(is.null(gene_list)|is.null(database_list)){
    stop('Source gene_list vector and database_list vector please')
  }
  #send a request with gene_list
  gene_list_collapsed=list(list=paste(gene_list,collapse='\n'))
  send=(httr::POST('http://amp.pharm.mssm.edu/Enrichr/addList',body=gene_list_collapsed)) #send a response
  user_list_id = httr::content(send,type='application/json')$userListId #retreive list id
  
  
  
  ###########extract GO data by libraries
  GO_enriched=list()
  
  for (j in database_list){
    
    raw_return=httr::GET(paste0('http://amp.pharm.mssm.edu/Enrichr/enrich?userListId=',user_list_id,'&backgroundType=',j))
    raw_return=httr::content(raw_return,type='application/json') #transform binary to list
    
    
    ################unlist obtained list
    genes=vector() #collapse genes into single field. 
	#some hardcoding stuff to collapse genes into single field or row 30-31 won't work and no parsing for obtained list. 
	# Gene names strickly placed in the [[x]][[y]][6] leaf of list, where [[x]]- name of the database and [[y]] GO entry (eg pathway).
    for (i in 1:length(raw_return[[j]])){
      genes=c(genes,paste(unlist(raw_return[[j]][[i]][6]),collapse=';')) 
      raw_return[[j]][[i]][6]=NULL
    }
    #### end of hardcoding
	
    #transform list to data.frame, rename columns and compute FDR since enlisted protocol does not gives it. WILL NOT WORK WO line 26-29 transformation!
    result_temp=data.frame(matrix(unlist(raw_return[[j]]), 
                                  nrow=length(raw_return[[j]]), byrow=TRUE)) 
    colnames(result_temp)= c('Rank', 'Term name', 'P-value', 'Z-score', 'Combined score', 'FDR', 'Old p-value', 'Old adjusted p-value')
    result_temp$genes=genes
#    result_temp$FDR=p.adjust(result_temp$`P-value`,method='fdr')
    
    if(output==TRUE){ #write result to working directory
      write.csv(result_temp,file=paste0(j,'_GO_enrichment.csv'))
    }
    
    GO_enriched[[j]]=result_temp
  } #End of database obtaining
  return(GO_enriched)
} #End of function
