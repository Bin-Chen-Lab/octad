#' @export
#' @import EnsDb.Hsapiens.v86
#' @import ensembldb


ensg_to_hgnc=function(data,select_surface=FALSE){
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= row.names(data), keytype = "GENEID", columns = c("SYMBOL","GENEID"))

data$GENEID=row.names(data)
data=merge(geneIDs1,data,by='GENEID')
data$GENEID=NULL
data=data[!duplicated(data$SYMBOL),]
row.names(data)=data$SYMBOL
data$SYMBOL=NULL
data=na.omit(data)
data$GenestableID=NULL
data$GENEID=NULL
if(select_surface==TRUE){
data=data[row.names(data)%in%basabsfinder::surface_expressed_genes$Gene,]
}
return(data)
}
