
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BioGA <a href="https://danymukesha.github.io/BioGA/"><img src="man/figures/logo.png" align="right" height="139" alt="BioGA website" /></a>

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
#> -- R CMD build ----------------------------------------------------------------------------------------------------------------
#>          checking for file 'C:\Users\dany.mukesha\AppData\Local\Temp\RtmpKsciOD\remotes35a879dd6114\danymukesha-BioGA-1917033/DESCRIPTION' ...  v  checking for file 'C:\Users\dany.mukesha\AppData\Local\Temp\RtmpKsciOD\remotes35a879dd6114\danymukesha-BioGA-1917033/DESCRIPTION' (369ms)
#>       -  preparing 'BioGA':
#>    checking DESCRIPTION meta-information ...     checking DESCRIPTION meta-information ...   v  checking DESCRIPTION meta-information
#> -  cleaning src
#>       -  checking for LF line-endings in source and make files and shell scripts
#>       -  checking for empty or unneeded directories
#>      Removed empty directory 'BioGA/vignettes'
#>    Omitted 'LazyData' from DESCRIPTION
#>       -  building 'BioGA_0.1.0.tar.gz'
#>      
#> 
#> Warning: package 'BioGA' is in use and will not be installed
```
