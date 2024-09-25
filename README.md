
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PCAone algorithms in R with RcppEigen\!

<!-- badges: start -->

<!-- badges: end -->

## Installation

``` r
# For the CRAN version
install.packages("pcaone")
# For the latest developing version
## devtools::install_github("Zilong-Li/PCAoneR")
```

## Example

This is a basic example which shows you how to use pcaone:

``` r
library(pcaone)
mat <- matrix(rnorm(100*5000), 100, 5000)
res <- pcaone(mat, k = 10)
str(res)
#> List of 3
#>  $ d: num [1:10] 80.1 79.3 78.8 78.5 78.4 ...
#>  $ u: num [1:100, 1:10] -0.282 -0.106 -0.0348 -0.0219 0.0414 ...
#>  $ v: num [1:5000, 1:10] -0.01971 0.00974 -0.02306 -0.00957 0.01311 ...
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
#>         expr       min        lq      mean    median        uq       max neval
#>          SVD 6.3386527 6.4493697 6.5878084 6.4936343 6.6752989 7.2448005    10
#>         rSVD 2.7598743 2.8006495 2.8523624 2.8390449 2.8630295 3.0286470    10
#>  pcaone.alg1 0.5111962 0.5174421 0.5360362 0.5257972 0.5529187 0.5814665    10
#>  pcaone.alg2 0.7594326 0.7668610 0.7872839 0.7833292 0.7878939 0.8441923    10
```

The above test is run on my MacBook Pro 2019 with processor 2.6 GHz
6-Core Intel Core i7. Note that the R is not linked to external BLAS or
MKL routine. To proper benchmark the performance with single core, we
can set the number of threads as one by `export OPENBLAS_NUM_THREADS=1
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1`.

## References

  - [Zilong Li, Jonas Meisner, Anders Albrechtsen (2023). Fast and
    accurate out-of-core PCA framework for large scale biobank
    data](https://genome.cshlp.org/content/33/9/1599)

## Todo

  - [ ] add `center` and `scale` method
