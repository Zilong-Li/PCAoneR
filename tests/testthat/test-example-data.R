library(pcaone)
#Set seed
set.seed(1234)

reconstruct <- function(s) s$u %*% diag(s$d) %*% t(s$v)

load(system.file("extdata", "popgen.rda", package="pcaone") )

## A <- popgen - rowMeans(popgen)
A <- sweep(popgen, 1, rowMeans(popgen), FUN = "-")
k <- 40

## A <- matrix(popgen, nrow = nrow(popgen))
d <- pcaone(A, k = k, opts = list("center" = TRUE, "scale" = TRUE, "byrow" = TRUE) )
## d <- pcaone(A, k = k, method = "winsvd", opts = list("center" = TRUE, "scale" = TRUE, "byrow" = TRUE) )

#*************************************************************************************
# Test: real tall matrix with population gentic data
#*************************************************************************************

system.time(s0 <- svd(A, nu = k, nv = k) )

## winSVD
system.time(s1 <- pcaone(A, k = k, method = "winsvd"))

max_d <- max(s0$d-s1$d)
Ae <- reconstruct(s1)
max_e <- max(abs(A-Ae))

testthat::test_that("Test: winSVD with nrow > ncol", {
  testthat::expect_true(max_d < 0.55)
  testthat::expect_true(max_e < 2.5)
})


system.time(s2 <- pcaone(A, k = k, method = "dashsvd"))

max_d <- max(s0$d-s2$d)
Ae <- reconstruct(s2)
max_e <- max(abs(A-Ae))

testthat::test_that("Test: dashSVD with nrow > ncol", {
  testthat::expect_true(max_d < 2.5)
  testthat::expect_true(max_e < 2.5)
})


#*************************************************************************************
# Test: real wide matrix with population gentic data
#*************************************************************************************

A <- t(A)
system.time(s0 <- svd(A, nu = k, nv = k) )

Ae <- s0$u %*% diag(s0$d[1:k]) %*% t(s0$v)

## winSVD
system.time(s1 <- pcaone(A, k = k, method = "winsvd"))

max_d <- max(s0$d-s1$d)
Ae <- reconstruct(s1)
max_e <- max(abs(A-Ae))

testthat::test_that("Test: winSVD with nrow > ncol", {
  testthat::expect_true(max_d < 0.55)
  testthat::expect_true(max_e < 2.5)
})


system.time(s2 <- pcaone(A, k = k, method = "dashsvd"))

max_d <- max(s0$d-s2$d)
Ae <- reconstruct(s2)
max_e <- max(abs(A-Ae))

testthat::test_that("Test: dashSVD with nrow > ncol", {
  testthat::expect_true(max_d < 2.5)
  testthat::expect_true(max_e < 2.5)
})


