library(pcaone)
data(popgen)
#Set seed
set.seed(1234)

A <- popgen - rowMeans(popgen)
k <- 40
system.time(s2 <- pcaone(A, k = k, method = "ssvd"))
system.time(s3 <- pcaone(A, k = k, method = "winsvd"))
system.time(s4 <- pcaone(A, k = k, method = "dashsvd"))
system.time(s5 <- PCAoneDashSVD(A, k = k, p = 7, q = 10, rand = 1))
system.time(s0 <- RSpectra::svds(A, k = k) )
s0$d-as.vector(s5$d)
