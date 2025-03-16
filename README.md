
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

  - `matrix` in base R
  - `dgeMatrix` in **Matrix** package, for general matrices
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
#>  $ d: num [1:10] 80.1 79.9 79.6 79 78.8 ...
#>  $ u: num [1:5000, 1:10] 0.00058 -0.02966 0.0031 0.0049 0.00569 ...
#>  $ v: num [1:100, 1:10] -0.3033 0.1388 -0.02 -0.1313 0.0323 ...
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
#>  27.474   3.009   1.279
system.time(s1 <- rsvd::rsvd(A, k = k, q = 4))  ## the number of epochs is two times of power iters, 4*2=8
#>    user  system elapsed 
#>   8.064  12.305   0.921
system.time(s2 <- pcaone(A, k = k, method = "ssvd", p = 7))   ## the number of epochs is 1 + p
#>    user  system elapsed 
#>   6.892   4.949   0.499
system.time(s3 <- pcaone(A, k = k, method = "winsvd", p = 7)) ## the number of epochs is 1 + p
#>    user  system elapsed 
#>   8.821   7.993   0.821
system.time(s4 <- pcaone(A, k = k, method = "dashsvd", p = 6))## the number of epochs is 2 + p
#>    user  system elapsed 
#>   5.290   1.041   0.295

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
#>  31.081  40.246   3.117
system.time(s2 <- pcaone(A, k = k, method = "ssvd", p = 20))
#>    user  system elapsed 
#>  16.221   8.828   1.064
system.time(s4 <- pcaone(A, k = k, method = "dashsvd", p = 18))
#>    user  system elapsed 
#>  13.323   2.598   0.697

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
print(timing, unit='s')
#> Unit: seconds
#>            expr       min        lq      mean    median        uq       max
#>        RSpectra 1.2188606 1.2569006 1.2934870 1.2716421 1.3021439 1.4671394
#>            rSVD 2.9463901 2.9751527 3.1790342 3.0345312 3.3626545 3.7373127
#>   pcaone.winsvd 0.7909559 0.8057960 0.8709076 0.8442994 0.9053499 1.1054273
#>     pcaone.ssvd 0.8526956 0.8768734 0.9126460 0.8968143 0.9404433 1.0432016
#>  pcaone.dashsvd 0.6797940 0.6874866 0.7363161 0.7036035 0.8150240 0.8522621
#>  neval
#>     10
#>     10
#>     10
#>     10
#>     10
```

## References

  - [Zilong Li, Jonas Meisner, Anders Albrechtsen (2023). Fast and
    accurate out-of-core PCA framework for large scale biobank
    data](https://genome.cshlp.org/content/33/9/1599)
  - [Feng et al. 2024. Algorithm 1043: Faster Randomized SVD with
    Dynamic Shifts](https://dl.acm.org/doi/10.1145/3660629)

## Todo

  - [ ] add `center` and `scale` method
