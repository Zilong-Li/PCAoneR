#' @title  Randomized Singular Value Decomposition Algorithm with Window-based Power Iterations from PCAone (Li et al 2022).
#
#' @description The Randomized Singular Value Decomposition (RSVD) computes the near-optimal low-rank approximation of a rectangular matrix
#' using a fast probablistic algorithm.
#
#' @details
#' The singular value decomposition (SVD) plays an important role in data analysis, and scientific computing.
#' Given a rectangular \eqn{(m,n)} matrix \eqn{A}, and a target rank \eqn{k << min(m,n)},
#' the SVD factors the input matrix \eqn{A} as
#'
#' \deqn{ A  =  U_{k} diag(d_{k}) V_{k}^\top }{ A = U diag(d) t(V)}
#'
#' The \eqn{k} left singular vectors are the columns of the
#' real or complex unitary matrix \eqn{U}. The \eqn{k} right singular vectors are the columns
#' of the real or complex unitary matrix \eqn{V}. The \eqn{k} dominant singular values are the
#' entries of \eqn{d}, and non-negative and real numbers.
#'
#' \eqn{q} is an oversampling parameter to improve the approximation.
#' A value of at least 10 is recommended, and \eqn{q=10} is set by default.
#'
#' The parameter \eqn{p} specifies the number of power (subspace) iterations
#' to reduce the approximation error. The power scheme is recommended,
#' especially when the singular values decay slowly. However, computing power iterations increases the
#' computational costs. Even though most RSVD implementations recommend \eqn{p=3} power iterations by default,
#' it's always sufficient to run only few power iterations where our window-based power iterations (\eqn{'alg2'})
#' come to play. We recommend using \eqn{windows=64} and \eqn{p>=7} for pcaone algorithm2. As it is designed for large dataset,
#' we recommend using \eqn{'alg2'} when \eqn{max(n,m) > 5000}.
#'
#' If \eqn{k > (min(n,m)/4)}, a deterministic partial or truncated \code{\link{svd}}
#' algorithm might be faster.
#'
#'
#' @param A       array_like; \cr
#'                a real/complex \eqn{(m, n)} input matrix (or data frame) to be decomposed.
#'
#' @param k       integer; \cr
#'                specifies the target rank of the low-rank decomposition. \eqn{k} should satisfy \eqn{k << min(m,n)}.
#'
#' @param p       integer, optional; \cr
#'                number of additional power iterations (by default \eqn{p=7}).
#'
#' @param s       integer, optional; \cr
#'                oversampling parameter (by default \eqn{s=10}).
#'
#' @param sdist   string \eqn{c( 'unif', 'normal')}, optional; \cr
#'                specifies the sampling distribution of the random test matrix: \cr
#'                		\eqn{'unif'} :  Uniform `[-1,1]`. \cr
#'                		\eqn{'normal'} (default) : Normal `~N(0,1)`. \cr
#'
#' @param method  string \eqn{c( 'winsvd', 'ssvd', 'winsvd')}, optional; \cr
#'                specifies the different variation of the randomized singular value decomposition : \cr
#'                		\eqn{'winsvd'} (default): window based RSVD in PCAone paper. \cr
#'                		\eqn{'ssvd'} : single-pass RSVD with power iterations in PCAone paper. \cr
#'                		\eqn{'dashsvd'} : RSVD with Dynamic Shifts by Feng et al. 2024. \cr
#'
#' @param batchs integer, optional; \cr
#'                the number of batchs for 'alg2' method. must be a power of 2 (by default \eqn{batchs=64}).
#'
#' @param shuffle logical, optional; \cr
#'                if shuffle the rows of input tall matrix or not for winsvd (by default \eqn{shuffle=TRUE}).
#'
#'@return \code{pcaone} returns a list containing the following three components:
#'\describe{
#'\item{d}{  array_like; \cr
#'           singular values; vector of length \eqn{(k)}.
#'}
#'
#'\item{u}{  array_like; \cr
#'           left singular vectors; \eqn{(m, k)} or \eqn{(m, nu)} dimensional array.
#'}
#'
#'\item{v}{  array_like; \cr
#'           right singular vectors; \eqn{(n, k)} or \eqn{(n, nv)} dimensional array. \cr
#'}
#'}
#'
#' @note The singular vectors are not unique and only defined up to sign.
#' If a left singular vector has its sign changed, changing the sign of the corresponding right vector
#' gives an equivalent decomposition.
#'
#'
#' @references
#' \itemize{
#'  \item Z. Li, J Meisner, A Albrechtsen. "Fast and accurate out-of-core PCA framework for large scale biobank data" (2023)
#'        \doi{10.1101/gr.277525.122}.
#' }
#'
#' @author Zilong Li \email{zilong.dk@gmail.com}
#'
#' @examples
#' library('pcaone')
#' data(popgen)
#' A <- popgen - rowMeans(popgen)
#' res <- pcaone(A, k = 10, p = 7, method = "winsvd")
#' str(res)
#' res <- pcaone(A, k = 10, p = 7, method = "ssvd")
#' str(res)
#' res <- pcaone(A, k = 10, p = 7, method = "dashsvd")
#' str(res)
#' @export
pcaone <- function(A, k=NULL, p=7, s=10, sdist="normal", method = "winsvd", batchs = 64, shuffle = TRUE) UseMethod("pcaone")

#' @export
pcaone.default <- function(A, k=NULL, p=7, s=10, sdist="normal", method = "winsvd", batchs = 64, shuffle = TRUE)
{
  rand <- switch(sdist,
                 normal = 1,
                 unif = 2,
                 stop("Selected sampling distribution is not supported!"))
  isPerm <- (shuffle && method == "winsvd")
  perm <- NULL
  isTall <- ifelse(nrow(A)>ncol(A), TRUE,  FALSE)
  
  if (isPerm) {
    n <- max(dim(A))
    perm <- sample(n)
    if(isTall) {
      A <- A[perm,]
    } else {
      A <- A[,perm]
    }
    
  }
  
  pcaoneObj <- switch(method,
                      ssvd = .Call(`_pcaone_PCAoneAlg1`, A, k, p, s, rand),
                      winsvd = .Call(`_pcaone_PCAoneAlg2`, A, k, p, s, rand, batchs),
                      dashsvd = dashSVD(A, k, p, s, rand),
                      stop("Method is not supported!"))
  pcaoneObj$d <- as.vector(pcaoneObj$d)
  pcaoneObj$u <- as.matrix(pcaoneObj$u)
  pcaoneObj$v <- as.matrix(pcaoneObj$v)

  if(!is.null(perm)) {
    original <- order(perm)
    if(isTall) {
      pcaoneObj$u <- pcaoneObj$u[original,]
    } else {
      pcaoneObj$v <- pcaoneObj$v[original,]
    }
  }

  class(pcaoneObj) <- "pcaone"
  return(pcaoneObj)
}

dashSVD <- function(A, k, p, s, rand) {
  if (nrow(A)>ncol(A)) {
    dashSVD_tall(A = A, k = k, p = p, s = s, rand = rand)
  } else {
    dashSVD_wide(A = A, k = k, p = p, s = s, rand = rand)
  }
}

#' @import stats
dashSVD_tall <- function(A, k, p = 3, s = 10, rand =  1) {
  M <- nrow(A)
  N <- ncol(A)
  ## L <- k + as.integer(ceiling(k / 2))
  L <- k + s
  if(rand == 1) {
    Omg <- matrix(stats::rnorm(L * M), M, L)
  } else {
    Omg <- matrix(stats::runif(L * M), M, L)
  }
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

#' @import stats
dashSVD_wide <- function(A, k, p = 3, s = 10, rand =  1) {
  M <- nrow(A)
  N <- ncol(A)
  ## L <- k + as.integer(ceiling(k / 2))
  L <- k + s
  if(rand == 1) {
    Omg <- matrix(stats::rnorm(L * N), N, L)
  } else {
    Omg <- matrix(stats::runif(L * N), N, L)
  }
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

# SVD via eigendecomposition 
eigSVD <- function(A, tol = 1e-10) {
  m <- nrow(A)
  n <- ncol(A)
  
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
  
  # Return results as list (thin SVD)
  list(U = U, S = sigma, V = V)
}
