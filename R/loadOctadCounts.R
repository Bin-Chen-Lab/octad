#' @export
#### runsRGES #######
loadOctadCounts <- function(sample_vector='',type='tpm',file=''){
  #  require(dplyr)
  #  require(RUVSeq)
  #  require(edgeR)
  transcripts = as.character(rhdf5::h5read(file, "meta/transcripts"))
  samples = as.character(rhdf5::h5read(file, "meta/samples"))
  if(type=='counts'){
    
    #print(paste('loading samples from whole OCTAD dataset log2 values for',length(sample_vector),'samples',sep=' '))
    print(paste('loading',length(transcripts),'log2 expression values for',length(sample_vector),'samples',sep=' '))
    
    exprData = rhdf5::h5read(file, "data/count", index=list(1:length(transcripts), which(samples %in% sample_vector)))
    rownames(exprData) = transcripts
    colnames(exprData) = samples[samples %in% sample_vector]
    rhdf5::H5close()
    return(exprData)
  }else if(type=='tpm'){
    
#    rownames(exprData) = transcripts
    print(paste('loading',length(transcripts),'TPM expression values for',length(sample_vector),'samples',sep=' '))
    
    exprData = rhdf5::h5read(file, "data/tpm", index=list(1:length(transcripts), which(samples %in% sample_vector)))
    rownames(exprData) = transcripts
    colnames(exprData) = samples[samples %in% sample_vector]
    
    rhdf5::H5close()
    return(exprData)
    #}else if(source=='octad.small'&type=='counts'){
    
    #print(paste('loading small OCTAD set containing only log2 expression values for 978 LINCS genes for',length(sample_vector),'samples',sep=' '))
    
    #exprData = rhdf5::h5read(file, "data/count", index=list(1:length(transcripts), which(samples %in% sample_vector)))
    #colnames(exprData) = samples[samples %in% sample_vector]
    #rownames(exprData) = transcripts
    #rhdf5::H5close()
    #return(exprData)
    #}else if(source=='octad.small'&type=='tpm'){
    
    #print(paste('loading small OCTAD set containing only TPM expression values for 978 LINCS genes for',length(sample_vector),'samples',sep=' '))
    
    #exprData = rhdf5::h5read(file, "data/tpm", index=list(1:length(transcripts), which(samples %in% sample_vector)))
    #colnames(exprData) = samples[samples %in% sample_vector]
    #rownames(exprData) = transcripts
    #rhdf5::H5close()
    #return(exprData)
  }else{
    print('Oops, some parameters are missing, check everything is ok and restart function')}
  rhdf5::H5close()
}