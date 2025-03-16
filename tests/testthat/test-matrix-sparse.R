library(pcaone)
library(Matrix)

m <- 5000
n <- 200
k <- 10

## Set up test matrices
set.seed(123)
x = matrix(rnorm(m * n), m)
x[sample(m * n, floor( m * n / 2))] = 0
# General matrices
gen = list(x,
           as(x, "dgeMatrix"),
           as(x, "dgCMatrix"),
           as(x, "dgRMatrix"))


# "true" values
s0 <- svd(x)


## Follow RSpectra
## Test whether the calculated (d, u, v) are consistent with svd()
## Return the largest residual
svd_resid <- function(res, svd0) {
  d_resid <- svd0$d[1:length(res$d)] - res$d
  u_resid = v_resid = 0
  if(!is.null(res$u))
    u_resid <- abs(svd0$u[, 1:ncol(res$u)]) - abs(res$u)
  if(!is.null(res$v))
    v_resid <- abs(svd0$v[, 1:ncol(res$v)]) - abs(res$v)
  mabs <- function(x) max(abs(x))
  maxerr <- max(mabs(d_resid), mabs(u_resid), mabs(v_resid))
  message(paste("max residual <", format(maxerr, digits = 5)))
  return(maxerr)
}


testthat::test_that("Test: winsvd for sparse matrix with nrow > ncol", {
  res <- sapply( lapply( gen, pcaone, k = k, method = "winsvd" ), svd_resid, svd0 = s0 )
  testthat::expect_identical(sum(as.vector(res) < 0.5), 4L)
})

testthat::test_that("Test: dashsvd for sparse matrix with nrow > ncol", {
  res <- sapply( lapply( gen, pcaone, k = k, method = "dashsvd" ), svd_resid, svd0 = s0 )
  testthat::expect_identical(sum(as.vector(res) < 1.0), 4L)
})

testthat::test_that("Test: ssvd for sparse matrix with nrow > ncol", {
  res <- sapply( lapply( gen, pcaone, k = k, method = "ssvd" ), svd_resid, svd0 = s0 )
  testthat::expect_identical(sum(as.vector(res) < 2.0), 4L)
})




