
library(pcaone)

#Set seed
set.seed(1234)

#*************************************************************************************
# Test: pcaone using real random test matrix
#*************************************************************************************
#Create real random test matrix of dimension m x n with target rank k
m = 500
n = 100
k = 10
p <- min(c(m, n))

L <- matrix(rnorm(m*m), m, m)
qrl <- qr(L)
U <- qr.Q(qrl)
L <- matrix(rnorm(n*n), n, n)
qrl <- qr(L)
V <- qr.Q(qrl)

# type2 matrix
## d <- (1:p)^(-2)
## testMat <- U %*% diag(d, m, n) %*% V


testMat <- matrix(runif(m*k), m, k)
testMat <- matrix(rnorm(m*k), m, k)
testMat <- testMat %*% t(testMat)
testMat <- testMat[,1:n]

#Deterministic SVD
svd_out <- svd(testMat)

#Randomized SVD k=n
rsvd_out <- pcaone(testMat, k = n, p = 1, method = "alg2")
testMat.re = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
testthat::test_that("Test 1: PCAoneLi k=n", {
              testthat::expect_equal(svd_out$d, rsvd_out$d)
              testthat::expect_equal(testMat, testMat.re)
          })

#Randomized SVD k=k
rsvd_out <- pcaone(testMat, k = k, windows = 2, p = 3, method = "alg2")
testMat.re = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
testthat::test_that("Test 2: PCAoneYu k=k", {
              testthat::expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
              testthat::expect_equal(testMat, testMat.re)
          })

#*************************************************************************************
# Test : testMat.T
#*************************************************************************************

testMat <- t(testMat)

#Deterministic SVD
svd_out <- svd(testMat)

#Randomized SVD k=n
rsvd_out <- pcaone(testMat, k = n, method = "alg2")
testMat.re = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
testthat::test_that("Test 1: PCAoneLi k=n", {
              testthat::expect_equal(svd_out$d, rsvd_out$d)
              testthat::expect_equal(testMat, testMat.re)
          })

#Randomized SVD k=k
rsvd_out <- pcaone(testMat, k = k, method = "alg2")
testMat.re = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
testthat::test_that("Test 2: PCAoneLi k=k", {
              testthat::expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
              testthat::expect_equal(testMat, testMat.re)
          })


