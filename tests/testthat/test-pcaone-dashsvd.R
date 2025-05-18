
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

#Randomized SVD k=k
rsvd_out <- pcaone(testMat, k = k, p = 1, s = 0, method = "dashsvd")
testMat.re <- rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
testthat::test_that("Test 2: dashSVD k=k, p=0, q=0", {
              testthat::expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
              testthat::expect_equal(testMat, testMat.re)
          })

#*************************************************************************************
# Test : testMat.T
#*************************************************************************************

testMat <- t(testMat)

#Deterministic SVD
svd_out <- svd(testMat)

#Randomized SVD k=k
rsvd_out <- pcaone(testMat, k = k, p = 2, s = 0, method = "dashsvd")
testMat.re = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
testthat::test_that("Test 2: dashSVD k=k, p=0, q=0", {
              testthat::expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
              testthat::expect_equal(testMat, testMat.re)
          })
