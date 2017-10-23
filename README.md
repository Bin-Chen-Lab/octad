Open Cancer TherApeutic Discovery (OCTAD)

When users specify a cancer or a list of clinical features, the system will retrieve relevant tumor sample profiles. A putative reference tissue will be inferred if not specified. A disease gene expression signature will be created, followed by gene set enrichment analysis. The signature will be later searched against the LINCS compound library, resulting in a ranked list of candidate compounds based on their reversal potency. Enriched targets/pathways of the drug hits will be performed. The disease signature will also be searched against the LINCS shRNA library, resulting in a list of candidate shRNAs based on their reversal potency. Enriched targets/pathways of the target hits will be performed. 


## to do list:
1. clean up the code. Run a few cases to ensure results are reproducible and meaniningful.
2. support query of clinical features (e.g., glioma samples with IDH mutant)
3. compound structure enrichment analysis, e.g.,  most interesting scaffolds of the top hits.
4. drug hits validation. One approach is to automatically collect drugs of the disease of interest from clinicaltrials.gov, and check their enrichment.
5. quality control. Some tumor samples that are impure or very different from cell lines, should be removed.
6. weight LINCS cell line. The RGES of each profile may be adjusted by the correlation of its cell line with the tumor samples of interest. We use gene expression to measure their correlation, we may add more.
7. support new RNA-Seq data (from EBI/GEO or user-generated).
8. parallel computing to speed up. Reduce the total running time to a few mins.

