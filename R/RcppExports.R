# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

svd_dense <- function(mat, k, p, s, batchs, params_pca) {
    .Call(`_pcaone_svd_dense`, mat, k, p, s, batchs, params_pca)
}

svd_sparse_col <- function(mat, k, p, s, batchs, params_pca) {
    .Call(`_pcaone_svd_sparse_col`, mat, k, p, s, batchs, params_pca)
}

svd_sparse_row <- function(mat, k, p, s, batchs, params_pca) {
    .Call(`_pcaone_svd_sparse_row`, mat, k, p, s, batchs, params_pca)
}

