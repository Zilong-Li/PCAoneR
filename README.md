
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Fast and accurate randomized singular value decomposition (RSVD)

## Introduction

There are 3 versions of RSVD implemented in this package, which are
ordered by their accuracy.

  - **winSVD: window based randomized singular value decomposition**
  - **dashSVD**: randomized singular value decomposition with dynamic
    shifted eigenvalues
  - **sSVD**: single pass randomized singular value decomposition with
    power iterations

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
mat <- matrix(rnorm(100*5000), 5000, 100)
res <- pcaone(mat, k = 10)
str(res)
#> List of 3
#>  $ d: num [1:10] 80.1 79.7 79.4 79 78.7 ...
#>  $ u: num [1:5000, 1:10] -0.00319 0.00769 0.01271 -0.00393 -0.00679 ...
#>  $ v: num [1:100, 1:10] 0.154527 0.000605 0.323732 -0.164047 0.052849 ...
#>  - attr(*, "class")= chr "pcaone"
```

## Accuracy

``` r
library(RSpectra) ## svds
library(rsvd)     ## regular rsvd
library(pcaone)
data(popgen)
A <- popgen - rowMeans(popgen) ## center
k <- 40
system.time(s0 <- RSpectra::svds(A, k = k) )
#>    user  system elapsed 
#>  29.057  11.899   1.765
system.time(s1 <- rsvd::rsvd(A, k = k, p = 7))
#>    user  system elapsed 
#>   5.647  12.359   0.795
system.time(s2 <- pcaone(A, k = k, method = "ssvd"))
#>    user  system elapsed 
#>  19.558  48.724   4.942
system.time(s3 <- pcaone(A, k = k, method = "winsvd"))
#>    user  system elapsed 
#>  40.248 117.970   9.466

par(mar = c(5, 5, 2, 2))
plot(s0$d-s1$d, ylim = c(0, 15), xlab = "PC index", ylab = "Error of singular valuse", cex = 2, cex.lab = 2)
points(s0$d-s2$d, col = "orange", cex = 2)
points(s0$d-s3$d, col = "red", cex = 2)
legend("top", legend = c("rSVD", "sSVD", "winSVD"), pch = 16,col = c("black", "orange", "red"), horiz = T, cex = 2, bty = "n" )
```

<img src="man/figures/README-acc-1.png" width="100%" />

## Benchmarking

Let’s see the performance of `pcaone` compared to the other rsvd
packages.

``` r
library(microbenchmark)
timing <- microbenchmark(
  'RSpectra' = svds(A,k = k),
  'rSVD' = rsvd(A, k=k, q = 7),
  'pcaone.rsvd' = pcaone(A, k=k, p = 7, method = "rsvd"),
  'pcaone.winsvd' = pcaone(A, k=k, p = 7),
  times=5)
print(timing, unit='s')
```

## References

  - [Zilong Li, Jonas Meisner, Anders Albrechtsen (2023). Fast and
    accurate out-of-core PCA framework for large scale biobank
    data](https://genome.cshlp.org/content/33/9/1599)

## Todo

  - [ ] add `center` and `scale` method
