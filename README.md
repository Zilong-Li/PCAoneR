
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Accurate and pass-efficient randomized SVD in R

## Introduction

This repo initially implements the algorithm so-called window-based
Randomized SVD in the [PCAone](https://github.com/Zilong-Li/PCAone)
paper. **Now the aim of the package is to implement state-of-the-art
Randomized SVD algorithms other than the basic
[rsvd](https://github.com/erichson/rSVD) for R community.**

Currently there are 2 versions of RSVD implemented in this package
ordered by their accuracy in below.

- **winSVD**: [window based randomized singular value
  decomposition](https://genome.cshlp.org/content/33/9/1599)
- **dashSVD**: [dynamic shifts based randomized singular value
  decomposition](https://dl.acm.org/doi/10.1145/3660629)

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
#>  $ d: num [1:10] 80.2 80 79.2 78.9 78.6 ...
#>  $ u: num [1:5000, 1:10] -0.01358 -0.01307 0.01905 0.00726 -0.01962 ...
#>  $ v: num [1:100, 1:10] 0.0516 -0.0284 0.0251 0.0852 -0.0507 ...
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
#>   4.781   0.012   4.797
system.time(s1 <- rsvd::rsvd(A, k = k, q = 4))  ## the number of epochs is two times of power iters, 4*2=8
#>    user  system elapsed 
#>   5.958   0.041   6.005
system.time(s3 <- pcaone(A, k = k, method = "winsvd", p = 7)) ## the number of epochs is 1 + p
#>    user  system elapsed 
#>   0.866   0.011   0.877
system.time(s4 <- pcaone(A, k = k, method = "dashsvd", p = 6))## the number of epochs is 2 + p
#>    user  system elapsed 
#>   0.668   0.025   0.694

par(mar = c(5, 5, 2, 1))
plot(s0$d-s1$d, ylim = c(0, 10), xlab = "PC index", ylab = "Error of singular values", cex = 1.5, cex.lab = 2)
points(s0$d-s3$d, col = "red", cex = 1.5)
points(s0$d-s4$d, col = "blue", cex = 1.5)
legend("top", legend = c("rSVD", "dashSVD", "winSVD"), pch = 16,col = c("black", "blue", "red"), horiz = T, cex = 1.2, bty = "n" )
```

<img src="man/figures/README-acc-1.png" width="100%" />

Now let’s see how many epochs we need for `rSVD`, `sSVD` and `dashSVD`
to reach the accuracy of `winSVD`.

``` r
system.time(s1 <- rsvd::rsvd(A, k = k, q = 20))  ## the number of epochs is 4*20=40
#>    user  system elapsed 
#>  24.976   0.158  25.165
system.time(s4 <- pcaone(A, k = k, method = "dashsvd", p = 18))
#>    user  system elapsed 
#>   1.745   0.054   1.799

par(mar = c(5, 5, 2, 1))
plot(s0$d-s1$d, ylim = c(0, 2), xlab = "PC index", ylab = "Error of singular values", cex = 1.5, cex.lab = 2)
points(s0$d-s3$d, col = "red", cex = 1.5)
points(s0$d-s4$d, col = "blue", cex = 1.5)
legend("top", legend = c("rSVD", "dashSVD", "winSVD"), pch = 16,col = c("black", "blue", "red"), horiz = T, cex = 1.2, bty = "n" )
```

<img src="man/figures/README-acc2-1.png" width="100%" /> \##
Benchmarking of speed

Let’s see the performance of `pcaone` compared to the other packages.

``` r
library(microbenchmark)
timing <- microbenchmark(
  'RSpectra' = svds(A,k = k),
  'rSVD' = rsvd(A, k=k, q = 20),
  'pcaone.winsvd' = pcaone(A, k=k, p = 7),
  'pcaone.dashsvd' = pcaone(A, k=k, p = 18, method = "dashsvd"),
  times=10)
#> Warning in microbenchmark(RSpectra = svds(A, k = k), rSVD = rsvd(A, k = k, :
#> less accurate nanosecond times to avoid potential integer overflows
print(timing, unit='s')
#> Unit: seconds
#>            expr        min         lq       mean     median         uq
#>        RSpectra  4.8384094  4.8984467  4.9434708  4.9301737  4.9701779
#>            rSVD 24.7174787 25.1742289 25.4543530 25.3989077 25.5714469
#>   pcaone.winsvd  0.8612599  0.8696085  0.8719048  0.8723306  0.8754653
#>  pcaone.dashsvd  1.7654236  1.7689123  1.7858041  1.7717733  1.7844873
#>         max neval
#>   5.1112929    10
#>  26.2576120    10
#>   0.8777382    10
#>   1.8411401    10
```

## References

- [Zilong Li, Jonas Meisner, Anders Albrechtsen (2023). Fast and
  accurate out-of-core PCA framework for large scale biobank
  data](https://genome.cshlp.org/content/33/9/1599)
- [Feng et al. 2024. Algorithm 1043: Faster Randomized SVD with Dynamic
  Shifts](https://dl.acm.org/doi/10.1145/3660629)
