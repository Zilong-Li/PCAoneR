
library(pcaone)

#Set seed
set.seed(1234)

#*************************************************************************************
# Test: pcaone using real random test matrix
#*************************************************************************************
#Create real random test matrix of dimension m x n with target rank k
m <- 50
n <- 30
k <- 10
testMat <- matrix(runif(m*k), m, k)
testMat <- testMat %*% t(testMat)
testMat <- testMat[,1:n]

#Deterministic SVD
svd_out <- svd(testMat)

#Randomized SVD k=n
rsvd_out <- pcaone(testMat, k = n, p = 0, method = "alg1")
testMat.re <- rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
testthat::test_that("Test 1: PCAoneYu k=n", {
              testthat::expect_equal(svd_out$d, rsvd_out$d)
              testthat::expect_equal(testMat, testMat.re)
          })

#Randomized SVD k=k
rsvd_out <- pcaone(testMat, k = k, p = 0, q = 0, method = "alg1")
testMat.re <- rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
testthat::test_that("Test 2: PCAoneYu k=k, p=0, q=0", {
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
rsvd_out <- pcaone(testMat, k = n, p = 3, q = 10, method = "alg1")
testMat.re = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
testthat::test_that("Test 1: PCAoneAlg1 k=n", {
              testthat::expect_equal(svd_out$d, rsvd_out$d)
              testthat::expect_equal(testMat, testMat.re)
          })

#Randomized SVD k=k
rsvd_out <- pcaone(testMat, k = k, p = 2, q = 0, method = "alg1")
testMat.re = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
testthat::test_that("Test 2: PCAoneAlg1 k=k, p=0, q=0", {
              testthat::expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
              testthat::expect_equal(testMat, testMat.re)
          })
