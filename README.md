
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

This is a basic example which shows you how to use pcaone:

``` r
library(pcaone)
mat <- matrix(rnorm(100*20000), 100, 20000)
res <- pcaone(mat, k = 10)
str(res)
#> List of 3
#>  $ d: num [1:10] 150 150 150 149 148 ...
#>  $ u: num [1:100, 1:10] 0.1192 -0.0258 0.1317 0.0323 0.1965 ...
#>  $ v: num [1:20000, 1:10] 0.008329 0.010267 -0.002753 0.000117 0.006507 ...
#>  - attr(*, "class")= chr "pcaone"
```
