
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Fast and accurate randomized singular value decomposition (RSVD)

## Introduction

There are 3 versions of RSVD implemented in this package, which are
ordered by their accuracy.

  - **winSVD**: [window based randomized singular value
    decomposition](https://genome.cshlp.org/content/33/9/1599)
  - **dashSVD**: [randomized singular value decomposition with dynamic
    shifted eigenvalues](https://dl.acm.org/doi/10.1145/3660629)
  - **sSVD**: [single pass randomized singular value decomposition with
    power iterations](https://genome.cshlp.org/content/33/9/1599)

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
#>  $ d: num [1:10] 81 80.1 79.5 79.1 78.6 ...
#>  $ u: num [1:5000, 1:10] -0.00402 -0.00015 -0.00481 0.01755 0.00365 ...
#>  $ v: num [1:100, 1:10] -0.1459 0.1601 0.0482 -0.0908 0.042 ...
#>  - attr(*, "class")= chr "pcaone"
```

## Benchmarking of accuracy

We define the accuracy as the error of singular values using results of
`RSpectra::svds` as truth. For all RSVD, let’s restrict the number of
epochs as 8 with only `winsvd` using 7, i.e. how many times the whole
matrix is read through if it can only be hold on disk.

``` r
library(RSpectra) ## svds
library(rsvd)     ## regular rsvd
library(pcaone)
data(popgen)
A <- popgen - rowMeans(popgen) ## center
k <- 40
system.time(s0 <- RSpectra::svds(A, k = k) )
#>    user  system elapsed 
#>  28.024  21.107   2.399
system.time(s1 <- rsvd::rsvd(A, k = k, q = 4))  ## the number of epochs is two times of power iters, 4*2=8
#>    user  system elapsed 
#>   9.027  20.213   1.881
system.time(s2 <- pcaone(A, k = k, method = "ssvd", p = 8))
#>    user  system elapsed 
#>   8.446   8.953   0.767
system.time(s3 <- pcaone(A, k = k, method = "winsvd", p = 7))
#>    user  system elapsed 
#>  10.177  10.926   1.024
system.time(s4 <- pcaone(A, k = k, method = "dashsvd", p = 6)) ## the number of epochs is 2 + power iters
#>    user  system elapsed 
#>  10.665  14.759   2.085

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
#>  36.119  69.747   5.707
system.time(s2 <- pcaone(A, k = k, method = "ssvd", p = 20))
#>    user  system elapsed 
#>  17.537  13.900   1.363
system.time(s4 <- pcaone(A, k = k, method = "dashsvd", p = 18))
#>    user  system elapsed 
#>  29.395  52.729   5.321

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
#>            expr       min       lq     mean   median       uq      max neval
#>        RSpectra 1.4525233 1.843104 2.050252 1.873966 2.306075 2.773527    10
#>            rSVD 4.5324825 4.658543 5.844413 5.714663 7.107667 7.930224    10
#>   pcaone.winsvd 1.0078065 1.079937 1.346716 1.217142 1.336223 2.296519    10
#>     pcaone.ssvd 0.9257167 1.059345 1.768647 1.867721 2.438595 2.569187    10
#>  pcaone.dashsvd 4.3526371 4.644177 4.820574 4.862795 5.048073 5.097065    10
```

## References

  - [Zilong Li, Jonas Meisner, Anders Albrechtsen (2023). Fast and
    accurate out-of-core PCA framework for large scale biobank
    data](https://genome.cshlp.org/content/33/9/1599)
  - [Feng et al. 2024. Algorithm 1043: Faster Randomized SVD with
    Dynamic Shifts](https://dl.acm.org/doi/10.1145/3660629)

## Todo

  - [ ] add `center` and `scale` method
