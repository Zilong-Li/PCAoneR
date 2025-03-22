
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Accurate and pass-efficient randomized SVD in R

## Introduction

This repo initially implements the algorithm so-called window-based
Randomized SVD in the [PCAone](https://github.com/Zilong-Li/PCAone)
paper. **Now the aim of the package is to implement state-of-the-art
Randomized SVD algorithms other than the basic
[rsvd](https://github.com/erichson/rSVD) for R community.**

Currently there are 3 versions of RSVD implemented in this package
ordered by their accuracy in below.

- **winSVD**: [window based randomized singular value
  decomposition](https://genome.cshlp.org/content/33/9/1599)
- **dashSVD**: [dynamic shifts based randomized singular value
  decomposition](https://dl.acm.org/doi/10.1145/3660629)
- **sSVD**: [single pass randomized singular value decomposition with
  power iterations](https://genome.cshlp.org/content/33/9/1599)

With surports for a number of matrix type including:

- `matrix` in base R for general dense matrices. or other types can be
  casted e.g. `dgeMatrix`
- `dgCMatrix` in **Matrix** package, for column major sparse matrices
- `dgRMatrix` in **Matrix** package, for row major sparse matrices

<!-- badges: start -->

<!-- badges: end -->

## Installation

``` r
# install.packages("pcaone") # For the CRAN version
remotes::install_github("Zilong-Li/PCAoneR") # For the latest developing version
```

## Example

This is a basic example which shows you how to use pcaone:

``` r
library(pcaone)
mat <- matrix(rnorm(100*5000), 5000, 100)
res <- pcaone(mat, k = 10)
str(res)
#> List of 3
#>  $ d: num [1:10] 80.2 79.8 79.3 78.9 78.7 ...
#>  $ u: num [1:5000, 1:10] -0.00311 0.00884 -0.01069 0.0206 -0.00773 ...
#>  $ v: num [1:100, 1:10] -0.02423 -0.08213 0.02153 0.00621 0.07479 ...
#>  - attr(*, "class")= chr "pcaone"
```

## Benchmarking of accuracy

We define the accuracy as the error of singular values using results of
`RSpectra::svds` as truth. For all RSVD, let’s restrict the number of
epochs as `8`, i.e. how many times the whole matrix is read through if
it can only be hold on disk.

``` r
library(RSpectra) ## svds
library(rsvd)     ## regular rsvd
library(pcaone)
load(system.file("extdata", "popgen.rda", package="pcaone") )
A <- popgen - rowMeans(popgen) ## center
k <- 40
system.time(s0 <- RSpectra::svds(A, k = k) )
#>    user  system elapsed 
#>   4.663   0.010   4.673
system.time(s1 <- rsvd::rsvd(A, k = k, q = 4))  ## the number of epochs is two times of power iters, 4*2=8
#>    user  system elapsed 
#>   5.892   0.033   6.195
system.time(s2 <- pcaone(A, k = k, method = "ssvd", p = 7))   ## the number of epochs is 1 + p
#>    user  system elapsed 
#>   0.697   0.020   0.717
system.time(s3 <- pcaone(A, k = k, method = "winsvd", p = 7)) ## the number of epochs is 1 + p
#>    user  system elapsed 
#>   0.823   0.012   0.834
system.time(s4 <- pcaone(A, k = k, method = "dashsvd", p = 6))## the number of epochs is 2 + p
#>    user  system elapsed 
#>   0.638   0.019   0.657

par(mar = c(5, 5, 2, 1))
plot(s0$d-s1$d, ylim = c(0, 10), xlab = "PC index", ylab = "Error of singular values", cex = 1.5, cex.lab = 2)
points(s0$d-s2$d, col = "orange", cex = 1.5)
points(s0$d-s3$d, col = "red", cex = 1.5)
points(s0$d-s4$d, col = "blue", cex = 1.5)
legend("top", legend = c("rSVD", "sSVD", "dashSVD", "winSVD"), pch = 16,col = c("black", "orange", "blue", "red"), horiz = T, cex = 1.2, bty = "n" )
```

<img src="man/figures/README-acc-1.png" width="100%" />

Now let’s see how many epochs we need for `rSVD`, `sSVD` and `dashSVD`
to reach the accuracy of `winSVD`.

``` r
system.time(s1 <- rsvd::rsvd(A, k = k, q = 20))  ## the number of epochs is 4*20=40
#>    user  system elapsed 
#>  24.600   0.115  24.717
system.time(s2 <- pcaone(A, k = k, method = "ssvd", p = 20))
#>    user  system elapsed 
#>   1.870   0.044   1.914
system.time(s4 <- pcaone(A, k = k, method = "dashsvd", p = 18))
#>    user  system elapsed 
#>   1.692   0.050   1.742

par(mar = c(5, 5, 2, 1))
plot(s0$d-s1$d, ylim = c(0, 2), xlab = "PC index", ylab = "Error of singular values", cex = 1.5, cex.lab = 2)
points(s0$d-s2$d, col = "orange", cex = 1.5)
points(s0$d-s3$d, col = "red", cex = 1.5)
points(s0$d-s4$d, col = "blue", cex = 1.5)
legend("top", legend = c("rSVD", "sSVD", "dashSVD", "winSVD"), pch = 16,col = c("black", "orange", "blue", "red"), horiz = T, cex = 1.2, bty = "n" )
```

<img src="man/figures/README-acc2-1.png" width="100%" />

## Benchmarking of speed

Let’s see the performance of `pcaone` compared to the other packages.

``` r
library(microbenchmark)
timing <- microbenchmark(
  'RSpectra' = svds(A,k = k),
  'rSVD' = rsvd(A, k=k, q = 20),
  'pcaone.winsvd' = pcaone(A, k=k, p = 7),
  'pcaone.ssvd' = pcaone(A, k=k, p = 20, method = "ssvd"),
  'pcaone.dashsvd' = pcaone(A, k=k, p = 18, method = "dashsvd"),
  times=10)
#> Warning in microbenchmark(RSpectra = svds(A, k = k), rSVD = rsvd(A, k = k, :
#> less accurate nanosecond times to avoid potential integer overflows
print(timing, unit='s')
#> Unit: seconds
#>            expr        min         lq       mean     median         uq
#>        RSpectra  4.6391060  4.6595379  4.7083525  4.6899738  4.7352504
#>            rSVD 24.3571732 24.3856144 24.4975993 24.4484750 24.4853403
#>   pcaone.winsvd  0.8603859  0.8650354  0.8687725  0.8688332  0.8714363
#>     pcaone.ssvd  1.9299811  1.9394407  1.9481555  1.9439759  1.9547296
#>  pcaone.dashsvd  1.7514591  1.7543524  1.7583468  1.7574516  1.7608106
#>         max neval
#>   4.9189169    10
#>  25.0868312    10
#>   0.8829495    10
#>   1.9822730    10
#>   1.7674187    10
```

## References

- [Zilong Li, Jonas Meisner, Anders Albrechtsen (2023). Fast and
  accurate out-of-core PCA framework for large scale biobank
  data](https://genome.cshlp.org/content/33/9/1599)
- [Feng et al. 2024. Algorithm 1043: Faster Randomized SVD with Dynamic
  Shifts](https://dl.acm.org/doi/10.1145/3660629)
