
library(devtools)
# export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
library(microbenchmark)
library(pcaone)
library(rsvd)

n <- 1000
m <- 5000
k <- 50
mat <- matrix(rnorm(n*m), n, m)

timing <- microbenchmark(
    'SVD' = svd(mat, nu=k, nv=k),
    'rSVD' = rsvd(mat, k=k, q = 7),
    'pcaone.alg1' = pcaone(mat, k=k, p = 7, method = "alg1"),
    'pcaone.alg2' = pcaone(mat, k=k, p = 7, method = "alg2"),
    times=10)
print(timing, unit = 's' )
