
# export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
library(microbenchmark)
library(pcaone)
library(rsvd)

n <- 2000
m <- 2000
k <- 50
mat <- matrix(rnorm(n*m), n, m)
res <- pcaone(mat, k = k, p = 3, q = 10, method = "li")

timing <- microbenchmark(
    'SVD' = svd(mat, nu = k, nv = k),
    'rSVD' = rsvd(mat, k = k, q = 3, p = 10),
    'PCAoneYu' = pcaone(mat, k = k, p = 3, q = 10, method = "yu"),
    'PCAoneLi' = pcaone(mat, k = k, p = 3, q = 10, method = "li"),
    times = 3)
print(timing, unit = 's' )
