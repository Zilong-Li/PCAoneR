
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PCAone algorithms in R with RcppEigen!

<!-- badges: start -->
<!-- badges: end -->

## Installation

You can install the development version of PCAoneR from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools") # if not installed
devtools::install_github("Zilong-Li/PCAoneR")
```

## Example

This is a basic example which shows you how to use PCAoneH:

``` r
library(pcaone)
mat <- matrix(rnorm(100*100), 100, 100)
res <- PCAoneY(mat, k = 10, l = 10, p = 3)
str(res)
#> List of 3
#>  $ d: num [1:10] 20.3 19 18.5 17.9 17.6 ...
#>  $ u: num [1:100, 1:10] -0.0556 0.224 0.0844 -0.0258 0.0499 ...
#>  $ v: num [1:100, 1:10] -0.037312 0.09289 0.080239 0.000949 -0.133278 ...
```
