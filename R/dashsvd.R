#' @title Dynamic shifted Randomized SVD 
#' @examples
#' mat <- matrix(rnorm(100*20000), 20000, 100)
#' res <- dashSVD(mat, k = 10)
#' str(res)
#' @export
dashSVD <- function(A, k, p = 10, s = 10) {
  if (nrow(A)>ncol(A)) {
    dashSVD_tall(A, k, p, s)    
  } else {
    dashSVD_wide(A, k, p, s)    
  }
}

dashSVD_tall <- function(A, k, p = 10, s = 10) {
  M <- nrow(A)
  N <- ncol(A)
  ## L <- k + as.integer(ceiling(k / 2))
  L <- k + s
  Omg <- matrix(stats::rnorm(L * M), M, L)
  Q <- t(A) %*% Omg
  e <- eigSVD(Q)
  Q <- e$U
  alpha <- 0.0
  for(i in 1:p) {
    ## message("power iteration ", i, ", alpha=",alpha)
    e <- eigSVD(t(A) %*% (A %*% Q)-alpha*Q)
    Q <- e$U
    if(e$S[L] > alpha) alpha <- (alpha + e$S[L]) / 2
  }
  e <- eigSVD(A %*% Q)
  U <- e$U
  S <- e$S
  V <- Q %*% e$V
  list(d = S[1:k], u = U[,1:k], v = V[,1:k])
}

dashSVD_wide <- function(A, k, p = 10, s = 10) {
  M <- nrow(A)
  N <- ncol(A)
  ## L <- k + as.integer(ceiling(k / 2))
  L <- k + s
  Omg <- matrix(stats::rnorm(L * N), N, L)
  Q <- A %*% Omg
  e <- eigSVD(Q)
  Q <- e$U
  alpha <- 0.0
  for(i in 1:p) {
    ## message("power iteration ", i, ", alpha=",alpha)
    e <- eigSVD(A %*% (t(A) %*% Q)-alpha*Q)
    Q <- e$U
    if(e$S[L] > alpha) alpha <- (alpha + e$S[L]) / 2
  }
  e <- eigSVD(t(A) %*% Q)
  V <- e$U
  S <- e$S
  U <- Q %*% e$V
  list(d = S[1:k], u = U[,1:k], v = V[,1:k])
}

# SVD via eigendecomposition for thin matrix 
eigSVD <- function(A, tol = 1e-10) {
  m <- nrow(A)
  n <- ncol(A)
  
  if (m >= n) {
    # Case 1: Compute via A^T A (smaller covariance matrix)
    C <- t(A) %*% A
    eig <- eigen(C, symmetric = TRUE)
    V <- eig$vectors
    sigma <- sqrt(pmax(eig$values, 0))  # Ensure non-negative singular values
    
    # Avoid division by near-zero values
    sigma_inv <- ifelse(sigma > tol, 1 / sigma, 0)
    U <- A %*% V %*% diag(sigma_inv)
    
    # Ensure U is orthonormal (QR decomposition)
    ## qrU <- qr(U)
    ## U <- qr.Q(qrU)
  } else {
    # Case 2: Compute via AA^T (smaller covariance matrix)
    C <- A %*% t(A)
    eig <- eigen(C, symmetric = TRUE)
    U <- eig$vectors
    sigma <- sqrt(pmax(eig$values, 0))  # Ensure non-negative singular values
    
    # Avoid division by near-zero values
    sigma_inv <- ifelse(sigma > tol, 1 / sigma, 0)
    V <- t(A) %*% U %*% diag(sigma_inv)
    
    # Ensure V is orthonormal (QR decomposition)
    ## qrV <- qr(V)
    ## V <- qr.Q(qrV)
  }
  
  # Return results as list (thin SVD)
  list(U = U, S = sigma, V = V)
}


