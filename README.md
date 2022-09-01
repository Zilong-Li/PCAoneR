
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
res <- pcaone(mat, k = 10, method = "li")
str(res)
#> List of 3
#>  $ d: num [1:10, 1] 150 150 149 149 149 ...
#>  $ u: num [1:100, 1:10] 0.1901 -0.0881 -0.1268 -0.0249 -0.0232 ...
#>  $ v: num [1:20000, 1:10] -0.000194 0.006939 -0.002693 -0.014784 0.002344 ...
#>  - attr(*, "class")= chr "pcaone"
```
