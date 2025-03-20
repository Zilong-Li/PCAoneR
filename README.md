
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
#>  $ d: num [1:10] 80 80 79.4 78.9 78.6 ...
#>  $ u: num [1:5000, 1:10] -0.0049 0.0192 -0.01689 0.00352 -0.00502 ...
#>  $ v: num [1:100, 1:10] -0.07187 -0.19122 -0.01924 -0.08518 0.00784 ...
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
#>  28.748   5.671   1.461
system.time(s1 <- rsvd::rsvd(A, k = k, q = 4))  ## the number of epochs is two times of power iters, 4*2=8
#>    user  system elapsed 
#>   8.630  13.607   1.017
system.time(s2 <- pcaone(A, k = k, method = "ssvd", p = 7))   ## the number of epochs is 1 + p
#>    user  system elapsed 
#>  17.099  34.516   3.563
system.time(s3 <- pcaone(A, k = k, method = "winsvd", p = 7)) ## the number of epochs is 1 + p
#>    user  system elapsed 
#>  37.849 108.382   7.642
system.time(s4 <- pcaone(A, k = k, method = "dashsvd", p = 6))## the number of epochs is 2 + p
#>    user  system elapsed 
#>   5.980   5.223   0.608

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
#>  32.674  45.492   3.436
system.time(s2 <- pcaone(A, k = k, method = "ssvd", p = 20))
#>    user  system elapsed 
#>  27.644  45.689   4.485
system.time(s4 <- pcaone(A, k = k, method = "dashsvd", p = 18))
#>    user  system elapsed 
#>  17.600  22.493   1.869

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
#>            expr      min       lq     mean   median       uq      max neval
#>        RSpectra 1.358759 1.380365 1.435317 1.399599 1.482455 1.623507    10
#>            rSVD 3.187601 3.323361 3.780321 3.444072 3.952823 6.216385    10
#>   pcaone.winsvd 7.563812 7.686590 8.090203 7.754446 7.763387 9.718475    10
#>     pcaone.ssvd 4.397840 4.474974 5.183663 4.576512 5.532823 8.190743    10
#>  pcaone.dashsvd 1.233052 1.280503 1.304921 1.309903 1.332356 1.354250    10
```

## References

  - [Zilong Li, Jonas Meisner, Anders Albrechtsen (2023). Fast and
    accurate out-of-core PCA framework for large scale biobank
    data](https://genome.cshlp.org/content/33/9/1599)
  - [Feng et al. 2024. Algorithm 1043: Faster Randomized SVD with
    Dynamic Shifts](https://dl.acm.org/doi/10.1145/3660629)

## Todo

  - [ ] add `center` and `scale` method
