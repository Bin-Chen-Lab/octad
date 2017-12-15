library(curl)

queryTCGA <- function(GENE, PROJECT){
#### GENE in ENSG format (e.g. ENSG00000138413)
#### PROJECT from TCGA (e.g. TCGA-LGG)
  
  ## retrieve all IDs with mutated GENE 
  all_gene<-curl_fetch_memory(paste0("https://api.gdc.cancer.gov/analysis/top_mutated_cases_by_gene?fields=submitter_id&pretty=true&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22genes.gene_id%22%2C%22value%22%3A%5B%22", GENE, "%22%5D%7D%7D%5D%7D&format=tsv&size=100000"))
  all_gene.df=read.table(text = rawToChar(all_gene$content), sep = '\t', header = TRUE)    
  all_gene.df$submitter_id=as.character(all_gene.df$submitter_id)
  
  ## retrieve all IDs with in PROJECT
  all_project<-curl_fetch_memory(paste0("https://api.gdc.cancer.gov/analysis/top_mutated_cases_by_gene?fields=submitter_id&pretty=true&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22",PROJECT,"%22%5D%7D%7D%5D%7D&format=tsv&size=100000"))
  all_project.df=read.table(text = rawToChar(all_project$content), sep = '\t', header = TRUE)    
  all_project.df$submitter_id=as.character(all_project.df$submitter_id)
  
  ## combine into grouped matrix
  combined=as.data.frame(unique(append(unique(all_project.df$submitter_id),unique(all_gene.df$submitter_id))))
  colnames(combined)[1]="submitter_id"
  
  # add in gene info
  combined$gene=0
  combined[which(combined$submitter_id %in% all_gene.df$submitter_id),"gene"]=1
  
  # add in project info
  combined$project=0
  combined[which(combined$submitter_id %in% all_project.df$submitter_id),"project"]=1  
  
  return(combined)
}





