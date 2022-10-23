
library(devtools)
# export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
library(microbenchmark)
library(pcaone)
library(rsvd)

n <- 2000
m <- 2000
k <- 50
mat <- matrix(rnorm(n*m), n, m)

res <- pcaone(mat, k = k, p = 3, q = 10, method = "yu", finder = 2)

timing <- microbenchmark(
    'SVD' = svd(mat, nu=k, nv=k),
    'rSVD' = rsvd(mat, k=k),
    'pcaone.alg1' = pcaone(mat, k=k, method = "alg1"),
    'pcaone.alg2' = pcaone(mat, k=k, method = "alg2"),
    times=10)
print(timing, unit = 's' )
