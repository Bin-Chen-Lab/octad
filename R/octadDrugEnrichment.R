#' @export
octadDrugEnrichment <- function(sRGES=NULL,target_type='chembl_targets',enrichFolder='enrichFolder'){
#  require(GSVA)
require(octad.db)
if(missing(sRGES)){
stop('sRGES input not found')
}
if(is.null(sRGES$sRGES)|is.null(sRGES$pert_iname)){
stop('Either sRGES or pert_iname collumn in Disease signature is missing')
}
options(warn=-1)
#  require(limma)
 if (!dir.exists(enrichFolder)) {
    dir.create(enrichFolder)
  }

for(target_type_selected in target_type){
cat(paste('Running enrichment for',target_type_selected,sep=' '),'\n')
  enrichFolder.n <- paste(enrichFolder,target_type_selected,sep='/')
  if (!dir.exists(enrichFolder.n)) {
    dir.create(enrichFolder.n)
  }
  
  #load random scores
#  load(paste0(dataFolder,"cmpd_sets_", target_type, ".RData"))
#  load(paste0(dataFolder,'random_gsea_score.RData'))
  
  #load required random scores from octad.db
  
  cmpd_sets = get(paste0("cmpd_sets_", target_type_selected), asNamespace('octad.db'))
  cmpdSets = cmpd_sets$cmpd.sets
  names(cmpdSets) = cmpd_sets$cmpd.set.names
  random_gsea_score=octad.db::random_gsea_score
  ############################
  
  
  drug_pred = sRGES
  
  rgess = matrix(-1*drug_pred$sRGES, ncol = 1)
  rownames(rgess) = drug_pred$pert_iname
  gsea_results = GSVA::gsva(rgess, cmpdSets, method = "ssgsea",  parallel.sz=8,ssgsea.norm = T,verbose=FALSE)
  
  gsea_results = merge(random_gsea_score[[target_type_selected]], gsea_results,by='row.names')
  row.names(gsea_results)=gsea_results$Row.names
  gsea_results$Row.names=NULL
  gsea_summary = data.frame(score = gsea_results[,ncol(random_gsea_score[[target_type_selected]])+1])
#calculating p.value
  gsea_p = apply(gsea_results, 1, function(x){
    sum(x[1:ncol(random_gsea_score[[target_type_selected]])] > x[ncol(random_gsea_score[[target_type_selected]])+1])/ncol(random_gsea_score[[target_type_selected]])
  })
  
  gsea_p = data.frame(target = names(gsea_p),score = gsea_summary, p = gsea_p, padj = p.adjust(gsea_p, method = 'fdr'))
  gsea_p = gsea_p[order(gsea_p$padj), ]
  # return(gsea_p)
  write.csv(gsea_p, paste0(enrichFolder.n, "/enriched_", target_type_selected, ".csv"),row.names = F)
  top.out.num = nrow(gsea_p[which(gsea_p$padj<=0.05),])
  if (top.out.num == 0) {
    top.out.num = 1
  }
  if(top.out.num > 50){
    top.out.num <- 50
  }
if(nrow(gsea_p)>0){
  for (i in 1:top.out.num) {
    top_target = as.character(gsea_p$target[i])
    sRGES$rank = rank(sRGES$sRGES)
    target_drugs_score = sRGES$rank[sRGES$pert_iname %in% cmpdSets[[top_target]]]
    if (length(target_drugs_score) < 3) {
      next
    }
    pdf(paste0(enrichFolder.n, "/top_enriched_", top_target, "_", target_type_selected, ".pdf"))
    limma::barcodeplot(sRGES$sRGES, target_drugs_score, main = top_target, xlab = "sRGES")
    dev.off()
  }
  if (target_type_selected=='ChemCluster'){
    clusternames <- as.character((gsea_p[which(gsea_p$padj<=0.05),])$target)
    if(length(clusternames) !=0){
    topclusterlist <- cmpdSets[clusternames]
    cat(sapply(topclusterlist, toString), file = paste0(enrichFolder.n,"misc.csv"), sep="\n")
    clusterdf <- read.csv2(paste0(enrichFolder.n,"misc.csv"), header=FALSE)
    clusterdf$cluster <- clusternames
    clusterdf$pval <- (gsea_p[which(gsea_p$padj<=0.05),])$padj
    colnames(clusterdf)[1] <- "drugs.in.cluster"
    write.csv(clusterdf,file=paste0(enrichFolder.n,'drugstructureclusters.csv'),row.names = F)
    }
  }
  cat(paste('Done for',target_type_selected,'for',nrow(gsea_p[which(gsea_p$padj<=0.05),]),'genes'),'\n')
  }else cat(paste('No signigicant enrichment found for',target_type_selected),'\n')
}
options(warn=0)
  }