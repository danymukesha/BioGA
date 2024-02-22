
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BioGA

<!-- badges: start -->

[![R-CMD-check](https://github.com/danymukesha/BioGA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/danymukesha/BioGA/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The `BioGA` package provides a set of functions for genetic algorithm
optimization tailored for analyzing high throughput genomic data. These
functions are implemented in C++ for improved speed and efficiency, with
an easy-to-use interface for use within R.

## Installation

You can install the package directly from GitHub using the `devtools`
package:

``` r
devtools::install_github("danymukesha/BioGA")
#> Downloading GitHub repo danymukesha/BioGA@HEAD
#> 
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#>          checking for file 'C:\Users\dany.mukesha\AppData\Local\Temp\RtmpqmJZpg\remotes921c790555c6\danymukesha-BioGA-fe9ca8f/DESCRIPTION' ...  ✔  checking for file 'C:\Users\dany.mukesha\AppData\Local\Temp\RtmpqmJZpg\remotes921c790555c6\danymukesha-BioGA-fe9ca8f/DESCRIPTION'
#>       ─  preparing 'BioGA':
#>    checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
#> ─  cleaning src
#>       ─  checking for LF line-endings in source and make files and shell scripts
#>       ─  checking for empty or unneeded directories
#>      Omitted 'LazyData' from DESCRIPTION
#>       ─  building 'BioGA_0.1.0.tar.gz'
#>      
#> 
#> Installing package into 'C:/Users/dany.mukesha/AppData/Local/R/win-library/4.3'
#> (as 'lib' is unspecified)
#> Warning in i.p(...): installation of package
#> 'C:/Users/DANY~1.MUK/AppData/Local/Temp/RtmpqmJZpg/file921c57703de/BioGA_0.1.0.tar.gz'
#> had non-zero exit status
```
