#specific key words, find their samples.
#source("https://bioconductor.org/biocLite.R")
#biocLite("TCGAbiolinks")
# project A list of valid project (it can be more than one) (see table below)
# data.category A valid project (see list with getProjectSummary(project))
# data.type A data type to filter the files to download
# sample.type A sample type to filter the files to download (See table below)
# workflow.type GDC workflow type
# barcode A list of barcodes to filter the files to download (can be partial barcodes)
# legacy Access legacy archive data (hg19 and hg18 data) instead of harmonized data? Default: FALSE
# platform Experimental data platform (HumanMethylation450, AgilentG4502A_07 etc). Used only for legacy repository
# file.type A string to filter files, based on its names. Used only for legacy repository

#http://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/query.html#useful_information

library(TCGAbiolinks)
