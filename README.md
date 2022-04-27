
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
res <- PCAoneH(mat, k = 10, l = 10, p = 3)
str(res)
#> List of 3
#>  $ d: num [1:10] 20.5 18.8 18.3 18.2 17.4 ...
#>  $ u: num [1:100, 1:10] -0.0988 0.0264 0.025 -0.0294 -0.1353 ...
#>  $ v: num [1:100, 1:10] 0.0176 0.0392 0.0515 0.1321 -0.089 ...
```
