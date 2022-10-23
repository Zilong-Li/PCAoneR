
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
#>  $ d: num [1:10] 149 148 148 148 147 ...
#>  $ u: num [1:100, 1:10] -0.0928 -0.1322 0.0149 -0.0666 0.049 ...
#>  $ v: num [1:20000, 1:10] 0.000918 0.002037 -0.016875 0.003456 0.018368 ...
#>  - attr(*, "class")= chr "pcaone"
```

## Benchmarking

Let’s see the performance of `pcaone` compared to the other rsvd
packages

``` r
library(microbenchmark)
library(pcaone)
library(rsvd)
data(tiger)
timing_svd <- microbenchmark(
    'SVD' = svd(tiger, nu=150, nv=150),
    'rSVD' = rsvd(tiger, k=150),
    'pcaone.alg1' = pcaone(tiger, k=150, method = "alg1"),
    'pcaone.alg2' = pcaone(tiger, k=150),
    times=10)
print(timing_svd, unit='s')
#> Unit: seconds
#>         expr       min        lq      mean    median        uq       max neval
#>          SVD 6.5842349 6.7716850 8.1404353 6.9975939 9.2145140 12.483671    10
#>         rSVD 2.1122792 2.1482307 2.3799463 2.2393917 2.3357059  3.188239    10
#>  pcaone.alg1 0.5377022 0.5515369 0.6812284 0.6225133 0.8445028  0.918856    10
#>  pcaone.alg2 0.7859225 0.7981073 0.9259427 0.8181026 0.8773590  1.412245    10
```
