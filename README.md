
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/danymukesha/BioGA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/danymukesha/BioGA/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

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
#> fs          (1.6.3   -> 1.6.4  ) [CRAN]
#> fastmap     (1.1.1   -> 1.2.0  ) [CRAN]
#> cachem      (1.0.8   -> 1.1.0  ) [CRAN]
#> xfun        (0.42    -> 0.44   ) [CRAN]
#> tinytex     (0.50    -> 0.51   ) [CRAN]
#> knitr       (1.45    -> 1.46   ) [CRAN]
#> htmltools   (0.5.7   -> 0.5.8.1) [CRAN]
#> bslib       (0.6.1   -> 0.7.0  ) [CRAN]
#> rmarkdown   (2.26    -> 2.27   ) [CRAN]
#> matrixStats (1.2.0   -> 1.3.0  ) [CRAN]
#> munsell     (0.5.0   -> 0.5.1  ) [CRAN]
#> farver      (2.1.1   -> 2.1.2  ) [CRAN]
#> BiocManager (1.30.22 -> 1.30.23) [CRAN]
#> bookdown    (0.38    -> 0.39   ) [CRAN]
#> gtable      (0.3.4   -> 0.3.5  ) [CRAN]
#> ggplot2     (3.5.0   -> 3.5.1  ) [CRAN]
#> Skipping 17 packages ahead of CRAN: BiocGenerics, graph, S4Arrays, IRanges, S4Vectors, MatrixGenerics, GenomeInfoDbData, zlibbioc, XVector, GenomeInfoDb, RBGL, Biobase, DelayedArray, GenomicRanges, BiocStyle, biocViews, SummarizedExperiment
#> Installing 16 packages: fs, fastmap, cachem, xfun, tinytex, knitr, htmltools, bslib, rmarkdown, matrixStats, munsell, farver, BiocManager, bookdown, gtable, ggplot2
#> Installing packages into 'C:/Users/dany.mukesha/AppData/Local/Temp/Rtmp63bptc/temp_libpath848868d23488'
#> (as 'lib' is unspecified)
#> Warning: unable to access index for repository https://bioconductor.org/packages/3.17/data/annotation/bin/windows/contrib/4.3:
#>   cannot open URL 'https://bioconductor.org/packages/3.17/data/annotation/bin/windows/contrib/4.3/PACKAGES'
#> Warning: unable to access index for repository https://bioconductor.org/packages/3.17/data/experiment/bin/windows/contrib/4.3:
#>   cannot open URL 'https://bioconductor.org/packages/3.17/data/experiment/bin/windows/contrib/4.3/PACKAGES'
#> Warning: unable to access index for repository https://bioconductor.org/packages/3.17/workflows/bin/windows/contrib/4.3:
#>   cannot open URL 'https://bioconductor.org/packages/3.17/workflows/bin/windows/contrib/4.3/PACKAGES'
#> package 'fs' successfully unpacked and MD5 sums checked
#> package 'fastmap' successfully unpacked and MD5 sums checked
#> package 'cachem' successfully unpacked and MD5 sums checked
#> package 'xfun' successfully unpacked and MD5 sums checked
#> package 'tinytex' successfully unpacked and MD5 sums checked
#> package 'knitr' successfully unpacked and MD5 sums checked
#> package 'htmltools' successfully unpacked and MD5 sums checked
#> package 'bslib' successfully unpacked and MD5 sums checked
#> package 'rmarkdown' successfully unpacked and MD5 sums checked
#> package 'matrixStats' successfully unpacked and MD5 sums checked
#> package 'munsell' successfully unpacked and MD5 sums checked
#> package 'farver' successfully unpacked and MD5 sums checked
#> package 'BiocManager' successfully unpacked and MD5 sums checked
#> package 'bookdown' successfully unpacked and MD5 sums checked
#> package 'gtable' successfully unpacked and MD5 sums checked
#> package 'ggplot2' successfully unpacked and MD5 sums checked
#> 
#> The downloaded binary packages are in
#>  C:\Users\dany.mukesha\AppData\Local\Temp\RtmpArcJMj\downloaded_packages
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#>          checking for file 'C:\Users\dany.mukesha\AppData\Local\Temp\RtmpArcJMj\remotes15e82c92423\danymukesha-BioGA-23ecb91/DESCRIPTION' ...  ✔  checking for file 'C:\Users\dany.mukesha\AppData\Local\Temp\RtmpArcJMj\remotes15e82c92423\danymukesha-BioGA-23ecb91/DESCRIPTION' (343ms)
#>       ─  preparing 'BioGA':
#>    checking DESCRIPTION meta-information ...     checking DESCRIPTION meta-information ...   ✔  checking DESCRIPTION meta-information
#> ─  cleaning src
#>       ─  checking for LF line-endings in source and make files and shell scripts
#>       ─  checking for empty or unneeded directories
#>      Omitted 'LazyData' from DESCRIPTION
#>       ─  building 'BioGA_0.99.5.tar.gz'
#>      
#> 
#> Installing package into 'C:/Users/dany.mukesha/AppData/Local/Temp/Rtmp63bptc/temp_libpath848868d23488'
#> (as 'lib' is unspecified)
```
