
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
mat <- matrix(rnorm(100*5000), 100, 5000)
res <- pcaone(mat, k = 10)
str(res)
#> List of 3
#>  $ d: num [1:10] 80.6 80.1 79.6 79.1 79 ...
#>  $ u: num [1:100, 1:10] -0.0565 -0.0404 -0.0268 -0.1161 0.0132 ...
#>  $ v: num [1:5000, 1:10] -0.001085 -0.01002 -0.001169 -0.000801 -0.015958 ...
#>  - attr(*, "class")= chr "pcaone"
```

## Benchmarking

Let’s see the performance of `pcaone` compared to the other rsvd
packages.

``` r
library(microbenchmark)
library(pcaone)
library(rsvd)
data(tiger)
timing <- microbenchmark(
    'SVD' = svd(tiger, nu=150, nv=150),
    'rSVD' = rsvd(tiger, k=150, q = 3),
    'pcaone.alg1' = pcaone(tiger, k=150, p = 3, method = "alg1"),
    'pcaone.alg2' = pcaone(tiger, k=150, p = 3, windows = 8),
    times=10)
print(timing, unit='s')
#> Unit: seconds
#>        expr       min        lq      mean    median        uq        max neval
#>         SVD 6.9045770 7.0373806 7.6021626 7.1162474 7.3132482 11.9759831    10
#>        rSVD 2.9189460 2.9601497 3.0412117 3.0452150 3.1174170  3.1540262    10
#> pcaone.alg1 0.5404496 0.5748938 0.6106437 0.5886755 0.6262602  0.8051892    10
#> pcaone.alg2 0.8177738 0.8211726 0.8621549 0.8587051 0.8740997  0.9599908    10
```

The above test is run on my MacBook Pro 2019 with processor 2.6 GHz
6-Core Intel Core i7. Note that the external BLAS or MKL routine is
disabled by
`export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1`.

## Todo

-   write `configure` to detect and use MKL
