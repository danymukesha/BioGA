
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
#> Using GitHub PAT from the git credential store.
#> Downloading GitHub repo danymukesha/BioGA@HEAD
#> Skipping 16 packages ahead of CRAN: BiocGenerics, graph, S4Arrays, IRanges, S4Vectors, MatrixGenerics, GenomeInfoDbData, zlibbioc, XVector, GenomeInfoDb, RBGL, Biobase, DelayedArray, GenomicRanges, biocViews, SummarizedExperiment
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#>          checking for file 'C:\Users\dany.mukesha\AppData\Local\Temp\RtmpeMTOmQ\remotes6e2c19f3527f\danymukesha-BioGA-39ae5be/DESCRIPTION' ...     checking for file 'C:\Users\dany.mukesha\AppData\Local\Temp\RtmpeMTOmQ\remotes6e2c19f3527f\danymukesha-BioGA-39ae5be/DESCRIPTION' ...   ✔  checking for file 'C:\Users\dany.mukesha\AppData\Local\Temp\RtmpeMTOmQ\remotes6e2c19f3527f\danymukesha-BioGA-39ae5be/DESCRIPTION' (387ms)
#>       ─  preparing 'BioGA':
#>    checking DESCRIPTION meta-information ...     checking DESCRIPTION meta-information ...   ✔  checking DESCRIPTION meta-information
#> ─  cleaning src
#>       ─  checking for LF line-endings in source and make files and shell scripts
#>       ─  checking for empty or unneeded directories
#>      Omitted 'LazyData' from DESCRIPTION
#>       ─  building 'BioGA_0.99.0.tar.gz'
#>      
#> 
#> Installing package into 'C:/Users/dany.mukesha/AppData/Local/Temp/RtmpQ9yEQT/temp_libpath43e05a565418'
#> (as 'lib' is unspecified)
```
