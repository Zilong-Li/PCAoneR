
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PCAoneR

<!-- badges: start -->
<!-- badges: end -->

PCAoneR brings PCAone algorithms in R by RcppEigen!

## Installation

You can install the development version of PCAoneR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Zilong-Li/PCAoneR")
```

## Example

This is a basic example which shows you how to use PCAoneH:

``` r
library(PCAoneR)
mat <- matrix(rnorm(100*100), 100, 100)
res <- PCAoneY(mat, k = 10, l = 10, p = 3)
str(res)
#> List of 3
#>  $ d: num [1:10] 19.4 18.8 18.3 18 17.9 ...
#>  $ u: num [1:100, 1:10] -0.0526 -0.0419 -0.0742 0.0993 -0.0708 ...
#>  $ v: num [1:100, 1:10] -0.0573 0.0566 0.2187 0.0994 0.0884 ...
```
