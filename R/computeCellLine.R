#' @export
#' @importFrom foreach foreach %do%
#' @importFrom rhdf5 h5read 


####### computeCellLine #######
computeCellLine <- function(case_id =case_id,
                                                        expSet=NULL,
                                                        LINCS_overlaps = TRUE,
                            source='octad.small',
                            file=NULL,
                                                        returnDF = FALSE){

    #STOPS
    if(missing(case_id)){
        stop('Case ids vector input not found')
    }




     #helper function by Ke
    pick.out.cell.line <- function(expr.of.samples,expr.of.cell.lines,marker.gene){
        marker.gene                     <- intersect(rownames(expr.of.samples),(marker.gene))
        marker.gene                     <- intersect(rownames(expr.of.cell.lines),(marker.gene))
        correlation.matrix        <- cor(expr.of.samples[marker.gene,],expr.of.cell.lines[marker.gene,],method='spearman')
        correlation.matrix[is.na(correlation.matrix)]=0
        cell.line.median.cor    <- apply(correlation.matrix,2,median) %>% sort(decreasing = TRUE)
        best.cell.line                <- names(cell.line.median.cor)[1]
        p.value.vec                     <- foreach::foreach(cell.line= setdiff(names(cell.line.median.cor),best.cell.line),.combine='c') %do% {
            v                                         <- correlation.matrix[,cell.line]
            p.value                             <- wilcox.test(correlation.matrix[,best.cell.line],v,alternative = 'greater',paired = TRUE)$p.value
        }
        names(p.value.vec) <- setdiff(names(cell.line.median.cor),best.cell.line)
        fdr.vec                        <- p.adjust(p.value.vec,method='fdr')
        list(cell.line.median.cor=cell.line.median.cor,best.cell.line=best.cell.line,compare.fdr.vec=fdr.vec,correlation.matrix = correlation.matrix )
    }
    #    load(paste0(dataFolder,'CCLE_OCTAD.RData'))
#    if(source=='octad.small'){
#        case_id=case_id[case_id%in%colnames(octad.db::octad.LINCS.counts)]
#        case_counts=octad.db::octad.LINCS.counts[,c(case_id)]
#    }else if (source=='octad.whole'){
if(source=='octad.whole'){
print(paste('loading whole octad expression data for',length(c(case_id)),'samples',sep=' '))

transcripts = as.character(rhdf5::h5read(file, "meta/transcripts"))
samples = as.character(rhdf5::h5read(file, "meta/samples"))
case_counts = rhdf5::h5read(file, "data/count", index=list(1:length(transcripts), which(samples %in% case_id)))
colnames(case_counts) = samples[samples %in% case_id]
rownames(case_counts) = transcripts
case_id = samples[samples %in% case_id]
#normal_counts = rhdf5::h5read(file, "data/count", index=list(1:length(transcripts), which(samples %in% control_id)))
#colnames(normal_counts) = samples[samples %in% control_id]
#rownames(normal_counts) = transcripts
#control_id = samples[samples %in% control_id]
H5close()

case_counts = cbind(case_counts)
#rm(normal_counts) # free up some memory
#rm(case_counts)
print(paste('computing correlation between cell lines and selected samples',sep=' '))
}else if(source=='octad.small'){
    print(paste('loading whole octad expression data for',length(c(case_id)),'samples',sep=' '))
    #H5close()
    case_counts=octad.db::octad.LINCS.counts[,c(case_id)]
    print(paste('computing correlation between cell lines and selected samples',sep=' '))
}else if(source!='octad'&missing(expSet)){
stop('Expression data not sourced, please, modify expSet option')

}else if(source != "octad.whole"|source != "octad.small"){
        case_counts=expSet[,case_id]
    }

    if(LINCS_overlaps == TRUE){

        CCLE.median                                 <- apply(octad.db::CCLE.overlaps,1,median)
    }else{
        CCLE.median = apply(octad.db::CCLE.log2.read.count.matrix,1,median)
    }
    CCLE.expressed.gene                 <- names(CCLE.median)[CCLE.median > 1]
    tmp                                                 <- octad.db::CCLE.log2.read.count.matrix[CCLE.expressed.gene,]
    tmp.rank                                        <- apply(tmp,2,rank)
    rank.mean                                     <- apply(tmp.rank,1,mean)
    rank.sd                                         <- apply(tmp.rank,1,sd)
    CCLE.rna.seq.marker.gene.1000                                 <- names(sort(rank.sd,decreasing =TRUE))[1:1000]
    TCGA.vs.CCLE.polyA.expression.correlation.result    <- pick.out.cell.line(case_counts, octad.db::CCLE.overlaps,CCLE.rna.seq.marker.gene.1000)
    correlation.dataframe <- TCGA.vs.CCLE.polyA.expression.correlation.result$cell.line.median.cor %>% as.data.frame()
    colnames(correlation.dataframe) <- "cor"
    cell.line.folder <- "CellLineEval/"
    if (!dir.exists(cell.line.folder)) {
        dir.create(cell.line.folder)
    }
    write.csv(correlation.dataframe, file = paste0(cell.line.folder,"CellLineCorrelations.csv"))

    topline <- data.frame(medcor = TCGA.vs.CCLE.polyA.expression.correlation.result$cell.line.median.cor) # could also do first
    # 3 of TCGA.vs.CCLE.polyA.expression.correlation.result$cell.line.median.cor


    if(returnDF == TRUE)
    {
        return(topline)
    }else{
        topline <- rownames(topline)[1]
        # topline$cellLine = row.names(topline)
        # topline = topline %>% select(cellLine,medcor)
        # load(paste0(dataFolder,'CCLE_OCTAD.RData'))
        # topline = left_join(topline,CCLE_demoDat %>% select(CCLE.shortname,Age,Gender,Race,Site.Primary,Histology,Hist.Subtype1),
        #    by=c('cellLine'='CCLE.shortname'))
        return(topline)
    }

}
