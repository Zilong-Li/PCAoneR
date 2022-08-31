
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
mat <- matrix(rnorm(100*200), 200, 100)
res <- pcaone(mat, k = 20, method = "yu")
str(res)
#> List of 3
#>  $ d: num [1:20, 1] 23.9 23.2 22.8 22.6 22.1 ...
#>  $ u: num [1:200, 1:20] 0.0678 0.0731 0.0793 0.0274 -0.0298 ...
#>  $ v: num [1:100, 1:20] -0.00622 0.0764 -0.0151 -0.04434 0.0464 ...
#>  - attr(*, "class")= chr "pcaone"
```
