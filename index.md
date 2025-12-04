# BioGA

*[BioGA](https://bioconductor.org/packages/3.20/BioGA)* package provides
a set of functions for genetic algorithm optimization adapted for
analyzing high throughput genomic data. These functions are implemented
in C++ for improved speed and efficiency, with an easy-to-use interface
for use within R.

## Installation

To install this package, start R (preferably version “4.4”) and enter:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install(pkgs = "BioGA", version = "devel", force = TRUE)
```

You can also install the package directly from GitHub using the
`devtools` package:

``` r
devtools::install_github("danymukesha/BioGA")
```
