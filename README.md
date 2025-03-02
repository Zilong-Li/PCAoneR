
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Fast and accurate randomized singular value decomposition (RSVD)

## Introduction

There are 3 versions of RSVD implemented in this package, which are
ordered by their accuracy.

  - **winSVD**: [window based randomized singular value
    decomposition](https://genome.cshlp.org/content/33/9/1599)
  - **dashSVD**: randomized singular value decomposition with dynamic
    shifted eigenvalues
  - **sSVD**: [single pass randomized singular value decomposition with
    power iterations](https://genome.cshlp.org/content/33/9/1599)

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
#>  $ d: num [1:10] 80.3 79.7 79.3 79.2 78.9 ...
#>  $ u: num [1:5000, 1:10] -0.00671 -0.01336 -0.0143 0.01195 0.00831 ...
#>  $ v: num [1:100, 1:10] 0.1466 0.1119 0.1024 0.2006 0.0736 ...
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
#>  27.916   1.570   1.237
system.time(s1 <- rsvd::rsvd(A, k = k, q = 4))
#>    user  system elapsed 
#>   9.673  18.263   1.645
system.time(s2 <- pcaone(A, k = k, method = "ssvd"))
#>    user  system elapsed 
#>  17.052  33.534   3.838
system.time(s3 <- pcaone(A, k = k, method = "winsvd"))
#>    user  system elapsed 
#>  40.115 120.547   9.508

par(mar = c(5, 5, 2, 2))
plot(s0$d-s1$d, ylim = c(0, 10), xlab = "PC index", ylab = "Error of singular values", cex = 2, cex.lab = 2)
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
