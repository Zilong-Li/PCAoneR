
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Accurate Randomized SVD with window-based power iteration

## Introduction

This repo initially implements the algorithm so-called window-based
Randomized SVD in the [PCAone](https://github.com/Zilong-Li/PCAone)
paper. **The aim of the package is to implement state-of-the-art
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

- `matrix` in base R for general dense matrices
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
#>  $ d: num [1:10] 80.4 79.6 79.4 78.8 78.8 ...
#>  $ u: num [1:5000, 1:10] -0.00662 0.01066 -0.02651 0.01647 0.00542 ...
#>  $ v: num [1:100, 1:10] 0.0128 -0.1413 0.033 -0.1198 -0.0105 ...
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
#>   4.719   0.012   4.733
system.time(s1 <- rsvd::rsvd(A, k = k, q = 4))  ## the number of epochs is two times of power iters, 4*2=8
#>    user  system elapsed 
#>   5.945   0.041   5.999
system.time(s2 <- pcaone(A, k = k, method = "ssvd", p = 7))   ## the number of epochs is 1 + p
#>    user  system elapsed 
#>   0.732   0.021   0.753
system.time(s3 <- pcaone(A, k = k, method = "winsvd", p = 7)) ## the number of epochs is 1 + p
#>    user  system elapsed 
#>   0.869   0.012   0.881
system.time(s4 <- pcaone(A, k = k, method = "dashsvd", p = 6))## the number of epochs is 2 + p
#>    user  system elapsed 
#>   0.670   0.025   0.695

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
#>  24.871   0.152  25.072
system.time(s2 <- pcaone(A, k = k, method = "ssvd", p = 20))
#>    user  system elapsed 
#>   1.919   0.045   1.965
system.time(s4 <- pcaone(A, k = k, method = "dashsvd", p = 18))
#>    user  system elapsed 
#>   1.738   0.059   1.797

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
#>            expr        min         lq       mean    median         uq
#>        RSpectra  4.7370844  4.7719100  4.8737781  4.913402  4.9283583
#>            rSVD 24.4685143 24.8870395 25.1024225 25.084026 25.2256136
#>   pcaone.winsvd  0.8668585  0.8703966  0.8771507  0.878141  0.8805611
#>     pcaone.ssvd  1.9454805  1.9495978  1.9992212  1.974393  2.0068966
#>  pcaone.dashsvd  1.7590869  1.7630143  1.7682542  1.767234  1.7730361
#>         max neval
#>   4.9584275    10
#>  25.6295077    10
#>   0.8913624    10
#>   2.1490443    10
#>   1.7807355    10
```

## References

- [Zilong Li, Jonas Meisner, Anders Albrechtsen (2023). Fast and
  accurate out-of-core PCA framework for large scale biobank
  data](https://genome.cshlp.org/content/33/9/1599)
- [Feng et al. 2024. Algorithm 1043: Faster Randomized SVD with Dynamic
  Shifts](https://dl.acm.org/doi/10.1145/3660629)
