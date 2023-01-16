# Web version:
http://octad.org/

# How to install
OCTAD and OCTAD.db is now available on the bioconductor:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("octad")
``` 

Latest vestion is available on the github:

``` 
library(devtools)
install_github("Bin-Chen-Lab/octad.db")
install_github("Bin-Chen-Lab/octad")
``` 

# Additional data
By default, octad package uses expression data for 978 genes from the LINCS dataset. However, it can influence the result and we advice using whole octad database. To obtatin whole results for DE, downloading of the additional OCTAD database [octad.counts.and.tpm.h5](https://experimenthub.bioconductor.org/fetch/7327) from the AWS link is required.


# Examples
Example workflow can be found in the package vignette:
```
vignette('octad')
``` 

# Contacts and citation
If you use our work, please cite the paper [OCTAD: an open workplace for virtually screening therapeutics targeting precise cancer patient groups using gene expression features, Nature Protocols](https://www.nature.com/articles/s41596-020-00430-z). Both OCTAD package and website was developed by [Bin Chen laboratory](http://binchenlab.org/).
Examples and questions can be addressed to Eugene Chekalin, PhD, chekali1@msu.edu or Bin Chen, PhD, PI, bin.chen@hc.msu.edu
