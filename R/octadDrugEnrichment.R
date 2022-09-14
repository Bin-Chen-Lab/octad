#' @export
#' @importFrom GSVA gsva
#' @importFrom limma barcodeplot
#' @import octad.db
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom S4Vectors mcols
#' @importFrom AnnotationHub query


octadDrugEnrichment <- function(sRGES = NULL, target_type = "chembl_targets", enrichFolder = "enrichFolder", outputFolder = NULL) {
  # require(GSVA)
  if (missing(sRGES)) {
    stop("sRGES input not found")
  }
  if (is.null(sRGES$sRGES) | is.null(sRGES$pert_iname)) {
    stop("Either sRGES or pert_iname collumn in Disease signature is missing")
  }
  options(warn = -1)

  if (is.null(outputFolder)) {
    outputFolder <- tempdir()
    message('outputFolder is NULL, writing output to tempdir()')
  }


  if (!dir.exists(enrichFolder)) {
    dir.create(file.path(outputFolder, enrichFolder))
  }
  #eh_dataframe <- as.data.frame(S4Vectors::mcols(AnnotationHub::query(.eh, "octad.db")))["title"]
  data('eh_dataframe',package='octad')
  random_gsea_score <- get_ExperimentHub_data("EH7275")
  for (target_type_selected in target_type) {
    message(paste("Running enrichment for", target_type_selected, sep = " "), "\n")
    enrichFolder.n <- file.path(outputFolder, enrichFolder, target_type_selected)
    if (!dir.exists(enrichFolder.n)) {
      dir.create(enrichFolder.n)
    }

    # load required random scores from octad.db
    eh_dataframe$object=row.names(eh_dataframe)
    cmpd_sets <- get_ExperimentHub_data((eh_dataframe[eh_dataframe$title == paste0("cmpd_sets_", target_type_selected),'object']))
    cmpdSets <- cmpd_sets$cmpd.sets
    names(cmpdSets) <- cmpd_sets$cmpd.set.names

    ############################


    drug_pred <- sRGES

    rgess <- matrix(-1 * drug_pred$sRGES, ncol = 1)

    if (is.null(dim(rgess))) {
      warning("rgess have zero rows, recompute it and try again")
      next
    } else {
      rownames(rgess) <- drug_pred$pert_iname
      rgess <- cbind(rgess, rgess) # PLUG FOR GSVA BUG, FIX AS THEY WOULD FIX THEIR CODE
      gsea_results <- GSVA::gsva(rgess, cmpdSets, method = "ssgsea", parallel.sz = 8, ssgsea.norm = TRUE, verbose = FALSE)
      gsea_results <- gsea_results[-1, ]
      gsea_results <- merge(random_gsea_score[[target_type_selected]], gsea_results, by = "row.names")

      if (is.null(dim(gsea_results))) {
        warning("gsea_results have zero rows, recompute it and try again")
        next
      } else {
        rownames(gsea_results) <- gsea_results$Row.names
        gsea_results$Row.names <- NULL
        gsea_summary <- data.frame(score = gsea_results[, ncol(random_gsea_score[[target_type_selected]]) + 1])
        # calculating p.value
        gsea_p <- apply(gsea_results, 1, function(x) {
          sum(x[seq_len(ncol(random_gsea_score[[target_type_selected]]))] > x[ncol(random_gsea_score[[target_type_selected]]) +
            1]) / ncol(random_gsea_score[[target_type_selected]])
        })
        gsea_p <- data.frame(target = names(gsea_p), score = gsea_summary, p = gsea_p, padj = p.adjust(gsea_p, method = "fdr"))
        gsea_p <- gsea_p[order(gsea_p$padj), ]
        write.csv(gsea_p, file.path(enrichFolder.n, paste0("enriched_", target_type_selected, ".csv")), row.names = FALSE)
        top.out.num <- nrow(gsea_p[which(gsea_p$padj <= 0.05), ])
        if (top.out.num == 0) {
          top.out.num <- 1
        }
        if (top.out.num > 50) {
          top.out.num <- 50
        }
        if (nrow(gsea_p) > 0) {
          for (i in seq_len(top.out.num)) {
            top_target <- as.character(gsea_p$target[i])
            sRGES$rank <- rank(sRGES$sRGES)
            target_drugs_score <- sRGES$rank[sRGES$pert_iname %in% cmpdSets[[top_target]]]
            if (length(target_drugs_score) < 3) {
              next
            }
            pdf(file.path(enrichFolder.n, paste0("/top_enriched_", top_target, "_", target_type_selected, ".pdf")))
            limma::barcodeplot(sRGES$sRGES, target_drugs_score, main = top_target, xlab = "sRGES")
            dev.off()
          }
          if (target_type_selected == "ChemCluster") {
            clusternames <- as.character((gsea_p[which(gsea_p$padj <= 0.05), ])$target)
            if (length(clusternames) != 0) {
              topclusterlist <- cmpdSets[clusternames]
              message(vapply(topclusterlist, toString), file = paste0(enrichFolder.n, "misc.csv"), sep = "\n")
              clusterdf <- read.csv2(paste0(enrichFolder.n, "misc.csv"), header = FALSE)
              clusterdf$cluster <- clusternames
              clusterdf$pval <- (gsea_p[which(gsea_p$padj <= 0.05), ])$padj
              colnames(clusterdf)[1] <- "drugs.in.cluster"
              write.csv(clusterdf, file = file.path(enrichFolder.n, "drugstructureclusters.csv"), row.names = FALSE)
            }
          }
          #msg <- paste("Done for", target_type_selected, "for", nrow(gsea_p[which(gsea_p$padj <= 0.05), ]), "genes")
          message("Done for", target_type_selected, "for", nrow(gsea_p[which(gsea_p$padj <= 0.05), ]), "genes")
        } else {
          #msg <- paste("No signigicant enrichment found for", target_type_selected)
          message("No signigicant enrichment found for", target_type_selected)
        }
      }
    }
  }
  options(warn = 0)
}
