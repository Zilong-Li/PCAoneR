
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Fast and accurate randomized singular value decomposition (RSVD)

## Introduction

There are 3 versions of RSVD implemented in this package, which are
ordered by their accuracy.

  - **winSVD**: [window based randomized singular value
    decomposition](https://genome.cshlp.org/content/33/9/1599)
  - **dashSVD**: [dynamic shifts based randomized singular
    value decomposition](https://dl.acm.org/doi/10.1145/3660629)
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
#>  $ d: num [1:10] 80.7 79.9 79.3 79.1 78.6 ...
#>  $ u: num [1:5000, 1:10] 0.00561 0.02096 0.02855 -0.00381 0.00558 ...
#>  $ v: num [1:100, 1:10] -0.02066 0.12137 0.11801 -0.11926 -0.00992 ...
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
data(popgen)
A <- popgen - rowMeans(popgen) ## center
k <- 40
system.time(s0 <- RSpectra::svds(A, k = k) )
#>    user  system elapsed 
#>  29.420  14.668   1.927
system.time(s1 <- rsvd::rsvd(A, k = k, q = 4))  ## the number of epochs is two times of power iters, 4*2=8
#>    user  system elapsed 
#>   8.722  17.691   1.317
system.time(s2 <- pcaone(A, k = k, method = "ssvd", p = 7))   ## the number of epochs is 1 + p
#>    user  system elapsed 
#>   7.122   6.239   0.775
system.time(s3 <- pcaone(A, k = k, method = "winsvd", p = 7)) ## the number of epochs is 1 + p
#>    user  system elapsed 
#>  10.122  11.114   1.034
system.time(s4 <- pcaone(A, k = k, method = "dashsvd", p = 6))## the number of epochs is 2 + p
#>    user  system elapsed 
#>   5.311   0.766   0.296

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
#>  34.010  56.014   4.426
system.time(s2 <- pcaone(A, k = k, method = "ssvd", p = 20))
#>    user  system elapsed 
#>  15.840   9.335   1.066
system.time(s4 <- pcaone(A, k = k, method = "dashsvd", p = 18))
#>    user  system elapsed 
#>  13.975   7.773   1.016

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
#>            expr       min        lq     mean   median       uq      max neval
#>        RSpectra 1.3123276 1.5685666 1.996660 1.888353 2.421958 2.975789    10
#>            rSVD 4.1148399 4.8318071 5.559569 5.266937 5.843289 8.291830    10
#>   pcaone.winsvd 0.8139188 0.8835514 1.247171 1.128622 1.401493 1.988634    10
#>     pcaone.ssvd 0.8720162 0.9468745 1.739007 1.736949 2.239249 3.085669    10
#>  pcaone.dashsvd 0.7368719 1.0958707 1.407987 1.403486 1.545074 2.092847    10
```

## References

  - [Zilong Li, Jonas Meisner, Anders Albrechtsen (2023). Fast and
    accurate out-of-core PCA framework for large scale biobank
    data](https://genome.cshlp.org/content/33/9/1599)
  - [Feng et al. 2024. Algorithm 1043: Faster Randomized SVD with
    Dynamic Shifts](https://dl.acm.org/doi/10.1145/3660629)

## Todo

  - [ ] add `center` and `scale` method
