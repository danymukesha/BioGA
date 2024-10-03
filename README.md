
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/danymukesha/BioGA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/danymukesha/BioGA/actions/workflows/R-CMD-check.yaml)
[![](https://img.shields.io/badge/devel%20version-0.99.6-blue.svg)](https://github.com/danymukesha/BioGA)
<!-- badges: end -->

# BioGA <a href="https://danymukesha.github.io/BioGA/"><img src="man/figures/logo.png" align="right" height="139" alt="BioGA website" /></a>

The `BioGA` package provides a set of functions for genetic algorithm
optimization tailored for analyzing high throughput genomic data. These
functions are implemented in C++ for improved speed and efficiency, with
an easy-to-use interface for use within R.

## Installation

To install this package, start R (version “4.4”) and enter folliw:

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
