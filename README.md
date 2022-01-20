# Open Cancer TherApeutic Discovery (OCTAD) database package

### Package overview
As the field of precision medicine progresses, we start to tailor treatments for cancer patients classified not only by their clinical, but also by their molecular features. The emerging cancer subtypes defined by these features require dedicated resources to assist the discovery of drug candidates for preclinical evaluation. Voluminous cancer patient gene expression profiles have been accumulated in public databases, enabling the creation of a cancer-specific expression signature. Meanwhile, large-scale gene expression profiles of chemical compounds became available in recent years. By matching the cancer-specific expression signature to compound-induced gene expression profiles from large drug libraries, researchers can prioritize small molecules that present high potency to reverse expression of signature genes for further experimental testing of their efficacy. This approach has proven to be an efficient and cost-effective way to identify efficacious drug candidates. However, the success of this approach requires multiscale procedures, imposing significant challenges to many labs. Therefore, we present OCTAD (http://octad.org): an open workspace for virtually screening compounds targeting precise cancer patient groups using gene expression features. We have included 19,127 patient tissue samples covering 50 cancer types, and expression profiles for 12,442 distinct compounds.  We will continue to include more tissue samples. We streamline all the procedures including deep-learning based reference tissue selection, disease gene expression signature creation, drug reversal potency scoring, and in silico validation. We release OCTAD as a web portal and a standalone R package to allow experimental and computational scientists to easily navigate the tool. The code and data can also be employed to answer various biological questions.

### How to install
To install the package run the following code:
```{r eval=FALSE} 
devtools::install_github('Bin-Chen-Lab/octad.db')
``` 
# Bioconductor installiation after the package will be released
To install the package run the following code:
```{r eval=FALSE} 
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("octad.db")
``` 

# Please note, this is a support package for the main package ```octad``` which can be obtained [there](https://github.com/Bin-Chen-Lab/octad)

  
  
###  Contacts and citation
  If you use our work, please cite the [OCTAD paper](https://www.nature.com/articles/s41596-020-00430-z). Both OCTAD package and website was developed by [Bin Chen laboratory](http://binchenlab.org/).
Examples and questions can be addressed to Evgenii Chekalin, PhD, chekali1@msu.edu or Bin Chen, PhD, PI, bin.chen@hc.msu.edu
