
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
#>  $ d: num [1:10] 81 80.3 79.5 79.5 79.3 ...
#>  $ u: num [1:5000, 1:10] 0.02376 0.01619 0.00688 0.01737 0.00316 ...
#>  $ v: num [1:100, 1:10] -0.0608 -0.0074 -0.0329 -0.0271 0.0213 ...
#>  - attr(*, "class")= chr "pcaone"
```

## Benchmarking of accuracy

We define the accuracy as the error of singular values using results of
`RSpectra::svds` as truth. For all RSVD, let’s restrict the number of
epochs as 8, i.e. how many times the whole matrix is read through if it
can only be hold on disk.

``` r
library(RSpectra) ## svds
library(rsvd)     ## regular rsvd
library(pcaone)
data(popgen)
A <- popgen - rowMeans(popgen) ## center
k <- 40
system.time(s0 <- RSpectra::svds(A, k = k) )
#>    user  system elapsed 
#>  28.247   1.710   1.258
system.time(s1 <- rsvd::rsvd(A, k = k, q = 4))  ## the number of epochs is two times of power iters, 4*2=8
#>    user  system elapsed 
#>   8.826  15.496   1.225
system.time(s2 <- pcaone(A, k = k, method = "ssvd", p = 8))
#>    user  system elapsed 
#>   6.660   5.916   0.954
system.time(s3 <- pcaone(A, k = k, method = "winsvd", p = 8))
#>    user  system elapsed 
#>  11.112   8.740   0.956
system.time(s4 <- pcaone(A, k = k, method = "dashsvd", p = 6)) ## the number of epochs is 2 + power iters
#>    user  system elapsed 
#>  10.654  12.905   1.792

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
#>  33.213  57.589   4.874
system.time(s2 <- pcaone(A, k = k, method = "ssvd", p = 20))
#>    user  system elapsed 
#>  15.453   5.560   0.881
system.time(s4 <- pcaone(A, k = k, method = "dashsvd", p = 20))
#>    user  system elapsed 
#>  30.136  37.707   4.358

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
  'pcaone.dashsvd' = pcaone(A, k=k, p = 20, method = "dashsvd"),
  times=10)
print(timing, unit='s')
#> Unit: seconds
#>            expr       min        lq     mean    median       uq      max neval
#>        RSpectra 1.2408488 1.4194163 1.766091 1.6250544 2.060196 2.908666    10
#>            rSVD 4.1430016 4.4327371 5.105631 4.9439714 5.785335 6.443811    10
#>   pcaone.winsvd 0.8368182 0.9249775 1.059625 0.9889866 1.076512 1.478980    10
#>     pcaone.ssvd 0.8181878 1.0357433 1.434357 1.0795606 1.906965 2.783159    10
#>  pcaone.dashsvd 4.7534786 5.0128576 5.295586 5.2492082 5.742555 5.798519    10
```

## References

  - [Zilong Li, Jonas Meisner, Anders Albrechtsen (2023). Fast and
    accurate out-of-core PCA framework for large scale biobank
    data](https://genome.cshlp.org/content/33/9/1599)
  - [Feng et al. 2024. Algorithm 1043: Faster Randomized SVD with
    Dynamic Shifts](https://dl.acm.org/doi/10.1145/3660629)

## Todo

  - [ ] add `center` and `scale` method
