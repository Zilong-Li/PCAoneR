
library(pcaone)

#Set seed
set.seed(1234)

#*************************************************************************************
# Test: pcaone using real random test matrix
#*************************************************************************************
#Create real random test matrix of dimension m x n with target rank k
m = 50
n = 30
k = 10
testMat <- matrix(runif(m*k), m, k)
testMat <- testMat %*% t(testMat)
testMat <- testMat[,1:n]

#Deterministic SVD
svd_out <- svd(testMat)

#Randomized SVD k=n
rsvd_out <- pcaone(testMat, k = n, method = "yu")
testMat.re = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)

testthat::test_that("Test 1: PCAoneYu k=n", {
              testthat::expect_equal(svd_out$d, rsvd_out$d)
              testthat::expect_equal(testMat, testMat.re)
          })

rsvd_out <- pcaone(testMat, k = n, p = 0, method = "li")
testMat.re = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)

testthat::test_that("Test 1: PCAoneLi k=n", {
              testthat::expect_equal(svd_out$d, rsvd_out$d)
              testthat::expect_equal(testMat, testMat.re)
          })
