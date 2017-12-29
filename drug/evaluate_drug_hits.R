#enrichment of known drugs

library("ROCR")
library("GSVA")
library("colorRamps")
dz_drugs <- read.csv(paste0("raw/validation/", dz, "_drugs.csv"), stringsAsFactors = F)

database <- "lincs"
sRGES = read.csv(paste0(dz, "/sRGES.csv"))
sRGES$drug <- 0
sRGES$drug[tolower(sRGES$pert_iname) %in% tolower(dz_drugs$Var1)] = 1

pred.roc=prediction(sRGES$sRGES,sRGES$drug)
performance(pred.roc,"auc")

sRGES <- sRGES[order(sRGES$sRGES), ]
sRGES[sRGES$drug == 1, ]

geneSets <- list(dz = sRGES[sRGES$drug == 1, "pert_iname"])

times <- 1000
ranks <- matrix(NA, nrow = nrow(sRGES), ncol = times)
ranks[,1] <- rank(sRGES$sRGES)
rownames(ranks) <- sRGES$pert_iname

for (i in 2:times){
  ranks[,i] <- rank(sample(1:length(sRGES$sRGES), length(sRGES$sRGES)))
}

gsea_results <- gsva(ranks, geneSets, method = "ssgsea")

print(paste("enriched p value", sum(gsea_results[1, -1] < gsea_results[1,1])/length(gsea_results[1, -1])))

enriched_p = sum(gsea_results[1, -1] < gsea_results[1,1])/length(gsea_results[1, -1])

par(mar=c(12, 1, 1, 1))

#adjust cmap score for better visualization
sRGES$score <- sapply(sRGES$sRGES, function(x){
  if (x < 0){
    -20^(abs(x))
  }else{
    20^(abs(x))
  }
})
pdf( paste0(dz, "/enrichment_", database, "_drugs.pdf", sep=""))
  par(mar=c(12, 4, 12, 4)) #bottom, left, top, and right
  z <- (matrix((sRGES$score)))
  n <- 200
  image(x=1:nrow(z), y = 1, z, col = green2red(n), xlab = "", ylab="", axes=FALSE )
  drugs_pos <- which(sRGES$drug == 1)
  for (drug_pos in drugs_pos){
    abline(v = drug_pos)
  }
  text(0.1, 0.1, "aaaaaa")
dev.off()

