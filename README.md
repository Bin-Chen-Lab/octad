# Web version:
http://octad.org/

# Bioconductor installation after the package will be released
To install the package run the following code:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("octad")
``` 

# How to Install
Before library installation install required octad.db, package with all required files for computation:
```
devtools::install_github('Bin-Chen-Lab/octad.db',build_vignettes = TRUE)
```
It takes a few minutes to install the package and verify files. Afterward, the pipeline will be ready to run. 
Install the package:
```
devtools::install_github('Bin-Chen-Lab/octad',build_vignettes = TRUE)
```



# Additional data
By default, octad package uses expression data for 978 genes from the LINCS dataset. However, it can influence the result and we advice using whole octad database. To obtatin whole results for DE, downloading of the additional OCTAD database [octad.counts.and.tpm.h5](https://chenlab-data-public.s3-us-west-2.amazonaws.com/octad/octad.counts.and.tpm.h5) from the AWS link is required.


# Examples
Example workflow can be found in the package vignette:
```
vignette('octad')
``` 

# Contacts and citation
If you use our work, please cite the paper [OCTAD: an open workplace for virtually screening therapeutics targeting precise cancer patient groups using gene expression features, Nature Protocols](https://www.nature.com/articles/s41596-020-00430-z). Both OCTAD package and website was developed by [Bin Chen laboratory](http://binchenlab.org/).
Examples and questions can be addressed to Eugene Chekalin, PhD, chekali1@msu.edu or Bin Chen, PhD, PI, bin.chen@hc.msu.edu
