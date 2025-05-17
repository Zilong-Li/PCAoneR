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
#'                the target rank for the low-rank decomposition. \eqn{k} should satisfy \eqn{k << min(m,n)}.
#'
#' @param p       integer, optional; \cr
#'                number of power iterations (by default \eqn{p=10}).
#'
#' @param s       integer, optional; \cr
#'                oversampling parameter (by default \eqn{s=10}).
#'
#' @param method  string \eqn{c( 'winsvd', 'dashsvd')}, optional; \cr
#'                specifies the different variation of the randomized singular value decomposition : \cr
#'                		\eqn{'auto'} (default): automatically choose the method based on size of input matrix. \cr
#'                		\eqn{'winsvd'} : window based RSVD in PCAone paper. \cr
#'                		\eqn{'dashsvd'} : dynamic shifts based RSVD with  by Feng et al. 2024. \cr
#'
#' @param B       integer, optional; \cr
#'                the number of batchs for 'winsvd' method. must be a power of 2 (by default \eqn{B=64}).
#'
#' @param shuffle logical, optional; \cr
#'                if shuffle the rows of input tall matrix or not for winsvd (by default \eqn{shuffle=TRUE}).
#'
#' @param opts    list, optional; \cr
#'                options related to PCA, e.g. center and scale. See details
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
#' load(system.file("extdata", "popgen.rda", package="pcaone") )
#' A <- popgen - rowMeans(popgen)
#' res <- pcaone(A, k = 40, method = "winsvd")
#' str(res)
#' res <- pcaone(A, k = 40, method = "dashsvd")
#' str(res)
#' @export
pcaone <- function(A, k=NULL, p=10, s=20, method = "auto", B = 64, shuffle = TRUE, opts = list()) UseMethod("pcaone")


#' @rdname pcaone
#' @export
pcaone.matrix <- function(A, k=NULL, p=10, s=20, method = "auto", B = 64, shuffle = TRUE, opts = list())
{
  ## A <- as.matrix(A)
  pcaopts <- check_pca_opts(A, opts)
  
  if(nrow(A)>1e4 || ncol(A)>1e4) {
    if(method == "auto")
      method <- "winsvd"
    else
      message("recommend using winsvd method for large matrix" )
  } else {
    if(method == "auto")
      method <- "dashsvd"
    else
      message("recommend using dashsvd method for small matrix" )
  }
  
  pcaopts$method <- switch(method,
                           winsvd = 1L,
                           dashsvd = 2L,
                           stop("Method is not supported!"))
  
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

  pcaoneObj <- .Call(`_pcaone_svd_dense`, A, k, p-1, s, B, pcaopts) ## 0-based
  pcaoneObj$d <- as.vector(pcaoneObj$d)
  pcaoneObj$u <- as.matrix(pcaoneObj$u)
  pcaoneObj$v <- as.matrix(pcaoneObj$v)
  if(pcaopts$dopca) {
    if(pcaopts$byrow)
      pcaoneObj$e <- pcaoneObj$d**2 / (ncol(A) - 1)
    else
      pcaoneObj$e <- pcaoneObj$d**2 / (nrow(A) - 1)
  }

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

#' @rdname pcaone
#' @export
pcaone.dgCMatrix <- function(A, k=NULL, p=10, s=20, method = "auto", B = 64, shuffle = TRUE, opts = list())
{
  pcaopts <- check_pca_opts(A, opts)

  if(nrow(A)>1e4 || ncol(A)>1e4) {
    if(method == "auto")
      method <- "winsvd"
    else
      message("recommend using winsvd method for large matrix" )
  } else {
    if(method == "auto")
      method <- "dashsvd"
    else
      message("recommend using dashsvd method for small matrix" )
  }

  pcaopts$method <- switch(method,
                           winsvd = 1L,
                           dashsvd = 2L,
                           stop("Method is not supported!"))
  
  pcaoneObj <- .Call(`_pcaone_svd_sparse_col`, A, k, p-1, s, B, pcaopts)
  pcaoneObj$d <- as.vector(pcaoneObj$d)
  pcaoneObj$u <- as.matrix(pcaoneObj$u)
  pcaoneObj$v <- as.matrix(pcaoneObj$v)
  if(pcaopts$dopca) {
    if(pcaopts$byrow)
      pcaoneObj$e <- pcaoneObj$d**2 / (ncol(A) - 1)
    else
      pcaoneObj$e <- pcaoneObj$d**2 / (nrow(A) - 1)
  }

  class(pcaoneObj) <- "pcaone"
  return(pcaoneObj)
}

#' @rdname pcaone
#' @export
pcaone.dgRMatrix <- function(A, k=NULL, p=7, s=10, method = "winsvd", B = 64, shuffle = TRUE, opts = list())
{
  pcaopts <- check_pca_opts(A, opts)

  if(nrow(A)>1e4 || ncol(A)>1e4) {
    if(method == "auto")
      method <- "winsvd"
    else
      message("recommend using winsvd method for large matrix" )
  } else {
    if(method == "auto")
      method <- "dashsvd"
    else
      message("recommend using dashsvd method for small matrix" )
  }

  pcaopts$method <- switch(method,
                           winsvd = 1L,
                           dashsvd = 2L,
                           stop("Method is not supported!"))
  
  pcaoneObj <- .Call(`_pcaone_svd_sparse_row`, A, k, p-1, s, B, pcaopts)
  pcaoneObj$d <- as.vector(pcaoneObj$d)
  pcaoneObj$u <- as.matrix(pcaoneObj$u)
  pcaoneObj$v <- as.matrix(pcaoneObj$v)
  if(pcaopts$dopca) {
    if(pcaopts$byrow)
      pcaoneObj$e <- pcaoneObj$d**2 / (ncol(A) - 1)
    else
      pcaoneObj$e <- pcaoneObj$d**2 / (nrow(A) - 1)
  }

  class(pcaoneObj) <- "pcaone"
  return(pcaoneObj)
}

check_pca_opts <- function(A, opts) {
  pcaopts <- list(rand = 1L,"dopca" = FALSE, "byrow" = FALSE, "center" = rep(0.0, 1), "scale" = rep(1.0, 1))
  if(is.null(opts$byrow)) opts$byrow <- FALSE
  pcaopts$byrow <- isTRUE(opts$byrow)
  if(is.null(opts$sdist)) opts$sdist <- "normal"
  pcaopts$rand <- switch(opts$sdist,
                         normal = 1L,
                         unif = 2L,
                         stop("Selected sampling distribution is not supported!"))

  if(isTRUE(opts$center)) {
    pcaopts$dopca <- TRUE
    ## center: either colMeans or rowMeans
    if(pcaopts$byrow) {
      pcaopts$center <-  rowMeans(A)
    } else {
      pcaopts$center <- colMeans(A)
    }
  } else if(is.numeric(opts$center)) {
    pcaopts$dopca <- TRUE
    n <- length(opts$center)
    if(pcaopts$byrow && n != nrow(A)) stop("opts$center must be TRUE/FALSE or a vector of length nrow(A) if byrow is TRUE")
    if(!pcaopts$byrow && n != ncol(A)) stop("opts$center must be TRUE/FALSE or a vector of length ncol(A) if byrow is FALSE")
    pcaopts$center <- opts$center
  } else {
    pcaopts$dopca <- FALSE
  }
  
  if(isTRUE(opts$scale)) {
    pcaopts$dopca <- TRUE
    if(pcaopts$byrow) {
      pcaopts$scale <- sqrt(rowSums(A**2) / (ncol(A)-1))
    } else {
      pcaopts$scale <- sqrt(colSums(A**2) / (nrow(A)-1))
    }
    pcaopts$scale[pcaopts$scale < 1e-9] <- 1.0
    ## instead of division, we multiply the inverse
    pcaopts$scale <- 1.0 / pcaopts$scale
  } else if(is.numeric(opts$scale)) {
    pcaopts$dopca <- TRUE
    n <- length(opts$scale)
    if(pcaopts$byrow && n != nrow(A)) stop("opts$center must be TRUE/FALSE or a vector of length nrow(A) if byrow is TRUE")
    if(!pcaopts$byrow && n != ncol(A)) stop("opts$center must be TRUE/FALSE or a vector of length ncol(A) if byrow is FALSE")
    pcaopts$scale <- opts$scale
  } else {
    pcaopts$dopca <- FALSE
  }

  if(pcaopts$dopca) {
    if(length(pcaopts$center) != length(pcaopts$center))
      stop("opts$center must has same length as scale" )
  }
  return(pcaopts)
}
