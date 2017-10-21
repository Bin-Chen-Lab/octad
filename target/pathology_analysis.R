#validate target
#progranosis

target = "CD143"
pathology = read.csv("raw/pathology.tsv", sep = "\t")

pathology[pathology$Gene.name %in% target,   ]
