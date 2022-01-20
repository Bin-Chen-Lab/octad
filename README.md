<<<<<<< HEAD
# Web version:
http://octad.org/

# How to Install
Before library installation install required Bioconductor and CRAN packages through this code:
```r
bioconductor_packages=c('edgeR','RUVSeq','DESeq2','limma','rhdf5','artMS')

#For R version 3.5> use BiocManager to install required bioconductor packages: 
if (length(setdiff(bioconductor_packages, rownames(installed.packages()))) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(setdiff(bioconductor_packages, rownames(installed.packages())))
}

packages=c('magrittr','dplyr','ggplot2','doParallel','foreach','lme4','Rfast','httr','data.table')
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
```
Next, install the octad.db, package with all required files for computation available via link  [octad.db](https://chenlab-data-public.s3.amazonaws.com/octad/octad.db_0.99.0.tar.gz%3Fdl%3D0)
```
devtools::install_github('Bin-Chen-Lab/octad.db',build_vignettes = TRUE)
```
It takes a few minutes to install the package and verify files. Afterward, the pipeline will be ready to run. 

```
install.packages("path%to%octad.db_0.99.0.tar.gz", repos = NULL, type="source")
```
Or without downloading the distributive:
```
install.packages("https://chenlab-data-public.s3-us-west-2.amazonaws.com/octad/octad.db_0.99.0.tar.gz",
                 method="libcurl",repos=NULL,type="source")
```
It takes a few minutes to install the package and verify files. Afterward, the pipeline will be ready to run. 
Finally, install the package:
```
devtools::install_github('Bin-Chen-Lab/octad',build_vignettes = TRUE)
```
=======
# Open Cancer TherApeutic Discovery (OCTAD) database package

### Package overview
As the field of precision medicine progresses, we start to tailor treatments for cancer patients classified not only by their clinical, but also by their molecular features. The emerging cancer subtypes defined by these features require dedicated resources to assist the discovery of drug candidates for preclinical evaluation. Voluminous cancer patient gene expression profiles have been accumulated in public databases, enabling the creation of a cancer-specific expression signature. Meanwhile, large-scale gene expression profiles of chemical compounds became available in recent years. By matching the cancer-specific expression signature to compound-induced gene expression profiles from large drug libraries, researchers can prioritize small molecules that present high potency to reverse expression of signature genes for further experimental testing of their efficacy. This approach has proven to be an efficient and cost-effective way to identify efficacious drug candidates. However, the success of this approach requires multiscale procedures, imposing significant challenges to many labs. Therefore, we present OCTAD (http://octad.org): an open workspace for virtually screening compounds targeting precise cancer patient groups using gene expression features. We have included 19,127 patient tissue samples covering 50 cancer types, and expression profiles for 12,442 distinct compounds.  We will continue to include more tissue samples. We streamline all the procedures including deep-learning based reference tissue selection, disease gene expression signature creation, drug reversal potency scoring, and in silico validation. We release OCTAD as a web portal and a standalone R package to allow experimental and computational scientists to easily navigate the tool. The code and data can also be employed to answer various biological questions.

### How to install
To install the package run the following code:
```{r eval=FALSE} 
devtools::install_github('Bin-Chen-Lab/octad.db')
``` 
>>>>>>> origin/master
# Bioconductor installiation after the package will be released
To install the package run the following code:
```{r eval=FALSE} 
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("octad.db")
``` 

<<<<<<< HEAD
# Additional data
By default, octad package uses expression data for 978 genes from the LINCS dataset. However, it can influence the result and we advice using whole octad database. To obtatin whole results for DE, downloading of the additional OCTAD database [octad.counts.and.tpm.h5](https://chenlab-data-public.s3-us-west-2.amazonaws.com/octad/octad.counts.and.tpm.h5) from the AWS link is required.

# Tutorial
The tutorial available via following [link](https://chenlab-data-public.s3-us-west-2.amazonaws.com/octad/octad_tutorial.pdf)

# Examples

The several examples listed in the file [octad_example.R](https://github.com/Bin-Chen-Lab/octad_desktop/blob/master/octad_example.R) :

<li>Example 1. liver hepatocellular carcinoma vs adjacent reference tissues;</li> 
<li>Example 2. breast cancer invasive carcinoma with PIK3 mutation vs reference tissues;</li> 
<li>Example 3. lung adenocarcinoma with amplified MYC gene vs reference tissues;</li> 
<li>Example 4. Primary breast cancer invasive carcinoma vs metastatic breast cancer invasive carcinoma using only 978 genes expression data for DE;</li> 
<li>Example 5. Compute sRGES score using GEO obtained dataset</li> 



# Contacts and citation
If you use our work, please cite the paper [OCTAD: an open workplace for virtually screening therapeutics targeting precise cancer patient groups using gene expression features, Nature Protocols](https://www.nature.com/articles/s41596-020-00430-z). Both OCTAD package and website was developed by [Bin Chen laboratory](http://binchenlab.org/).
Examples and questions can be addressed to Eugene Chekalin, PhD, chekali1@msu.edu or Bin Chen, PhD, PI, bin.chen@hc.msu.edu
=======
# Please note, this is a support package for the main package ```octad``` which can be obtained [there](https://github.com/Bin-Chen-Lab/octad)

  
  
###  Contacts and citation
  If you use our work, please cite the [OCTAD paper](https://www.nature.com/articles/s41596-020-00430-z). Both OCTAD package and website was developed by [Bin Chen laboratory](http://binchenlab.org/).
Examples and questions can be addressed to Evgenii Chekalin, PhD, chekali1@msu.edu or Bin Chen, PhD, PI, bin.chen@hc.msu.edu
>>>>>>> origin/master
