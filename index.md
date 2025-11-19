# `CroCoNet` Cross-species Comparison of Networks

CroCoNet is a framework to quantitatively compare gene regulatory
networks across species and identify conserved and diverged modules. It
hinges on contrasting module variability within and across species in
order to distinguish truw evolutionary divergence from confounding
factors such as diversity across individuals, environmental differences
and technical noise.

![Pipeline](reference/figures/pipeline.png)  
^(The main steps of the CroCoNet workflow)

## ⬇️ Installation

For the installation, the R package `devtools` is needed.

``` r
install.packages("devtools")
library(devtools)
```

Once `devtools` is available, you can install the development version of
`CroCoNet` and all its dependencies from GitHub with:

``` r
devtools::install_github("Hellmann-Lab/CroCoNet")
```

## 📖 User guide

For a step-by-step guide and detailed explanations, please check out the
vignette on the analysis of an example scRNA-seq dataset:

``` r
browseVignettes("CroCoNet")
```

You can access this vignette, along with the documentation of all
functions, on the [CroCoNet
website](https://hellmann-lab.github.io/CroCoNet/hellmann-lab.github.io/CroCoNet/)
as well.

## 📜 Citation

The CroCoNet framework and its biological applications are described in
detail in our
[preprint](https://www.biorxiv.org/content/10.1101/2025.11.18.689002v1)
on bioRxiv. Please cite this manuscript as follows:

``` r
citation("CroCoNet")
```
