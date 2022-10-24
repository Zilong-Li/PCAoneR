
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
#>  $ d: num [1:10] 149 149 148 148 147 ...
#>  $ u: num [1:100, 1:10] 0.0314 -0.057 0.0564 -0.2068 0.1022 ...
#>  $ v: num [1:20000, 1:10] 0.008411 -0.000617 -0.00383 0.000789 -0.003622 ...
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
timing_svd <- microbenchmark(
    'SVD' = svd(tiger, nu=150, nv=150),
    'rSVD' = rsvd(tiger, k=150, q = 3),
    'pcaone.alg1' = pcaone(tiger, k=150, method = "alg1"),
    'pcaone.alg2' = pcaone(tiger, k=150, method = "alg2"),
    times=10)
print(timing_svd, unit='s')
#> Unit: seconds
#>         expr       min        lq      mean    median        uq        max neval
#>          SVD 6.5416716 6.5732372 7.2773343 6.7346177 6.8681427 11.7746731    10
#>         rSVD 2.8050177 2.9114025 2.9568331 2.9442578 2.9986419  3.1445448    10
#>  pcaone.alg1 0.5164194 0.5219948 0.5395647 0.5249782 0.5373050  0.6395444    10
#>  pcaone.alg2 0.7732617 0.7956074 0.8242164 0.8260486 0.8518961  0.8903296    10
```

The above test is run on my MacBook Pro 2019 with processor 2.6 GHz
6-Core Intel Core i7. Note that the external BLAS or MKL routine is
disabled by
`export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1`.

## Todo

-   write `configure` to detect and use MKL
