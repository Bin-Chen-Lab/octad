library(curl)
library(estimate)

id_mapping <-function(id, input_format = "hgnc_symbol", output_format = "ensembl_gene_id"){
  #ensembl_gene_id ensembl_transcript_id hgnc_symbol entrezgene
  library(biomaRt)
  
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  id_mapped <- getBM(attributes=c(input_format, output_format), filters =
                       input_format, values =id, mart = ensembl)
  return(id_mapped[1,2]) #only return the first hit
}

id_mapping_gene_ensembl <-function(id){
  #ensembl_gene_id ensembl_transcript_id hgnc_symbol entrezgene
  mapping = read.csv("raw/gencode.v23.annotation.gene.probeMap.csv", stringsAsFactors = F)
  return(mapping$ensembl[mapping$gene == id][1])
}

queryGDC <- function(GENE, PROJECT){
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


estimatePurity  <- function(expr_matrix){
  #expr_matrix: with gene symbols as row names and samples as colnames
  #return purity score
  samples = data.frame(NAME = rownames(expr_matrix), Description = NA,  expr_matrix)
  write("", file = "temp_samples.gct")
  write(paste(nrow(samples), "\t", ncol(samples) -2), file = "temp_samples.gct", append = T)
  write("", file = "temp_samples.gct", append = T)
  write.table(samples,  file = "temp_samples.gct", append = T, row.names=F, col.names=T, sep="\t", quote=F)
  
  in.file <- ("temp_samples.gct")
  out.file <- "temp_samples_output.gct"
  estimateScore(in.file, out.file)
  
  estimateScore = read.delim(out.file, sep = "\t", skip = 3)
  
  return(as.numeric(estimateScore[3, -c(1,2)]))
}



