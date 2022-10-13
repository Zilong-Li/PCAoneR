#' @title  Randomized Singular Value Decomposition Algorithms in PCAone (pcaone).
#
#' @description The randomized SVD computes the near-optimal low-rank approximation of a rectangular matrix
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
#' \eqn{p} is an oversampling parameter to improve the approximation.
#' A value of at least 10 is recommended, and \eqn{p=10} is set by default.
#'
#' The parameter \eqn{q} specifies the number of power (subspace) iterations
#' to reduce the approximation error. The power scheme is recommended,
#' if the singular values decay slowly. In practice, 2 or 3 iterations
#' achieve good results, however, computing power iterations increases the
#' computational costs. The power scheme is set to \eqn{q=2} by default.
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
#'                number of additional power iterations (by default \eqn{p=3}).
#'
#' @param q       integer, optional; \cr
#'                oversampling parameter (by default \eqn{q=10}).
#'
#' @param windows integer; \cr
#'                the number of windows for 'alg2' method. must be a power of 2 (by default \eqn{windows=2^{p+2}}).
#'
#' @param sdist   string \eqn{c( 'unif', 'normal')}, optional; \cr
#'                specifies the sampling distribution of the random test matrix: \cr
#'                		\eqn{'unif'} :  Uniform `[-1,1]`. \cr
#'                		\eqn{'normal'} (default) : Normal `~N(0,1)`. \cr
#'
#' @param method  string \eqn{c( 'alg1', 'agl2')}, optional; \cr
#'                specifies the different variation of the randomized singular value decomposition : \cr
#'                		\eqn{'alg2'} (default): PCAone: fast and accurate out-of-core PCA framework for large scale biobank data (2022). \cr
#'                		\eqn{'alg1'} : Single-Pass PCA of Large High-Dimensional Data (2017). \cr
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
#' @note The singular vectors are not unique and only defined up to sign
#' (a constant of modulus one in the complex case). If a left singular vector
#' has its sign changed, changing the sign of the corresponding right vector
#' gives an equivalent decomposition.
#'
#'
#' @references
#' \itemize{
#'  \item [1] Z. Li, J Meisner, A Albrechtsen. "PCAone: fast and accurate out-of-core PCA framework for large scale biobank data" (2022)
#'        \url{https://doi.org/10.1101/2022.05.25.493261}.
#' }
#'
#' @author Zilong Li \email{zilong.dk@gmail.com}
#'
#' @examples
#'library('pcaone')
#' mat <- matrix(rnorm(100*20000), 100, 20000)
#' res <- pcaone(mat, k = 10, p = 5, windows = 32, method = "li")
#' str(res)
#' res <- pcaone(mat, k = 10, p = 5, method = "yu")
#' str(res)
#' @export
pcaone <- function(A, k=NULL, windows = NULL, p=3, q=10, sdist="normal", method = "alg2", finder = 1) UseMethod("pcaone")

#' @export
pcaone.default <- function(A, k=NULL, windows = NULL, p=3, q=10, sdist="normal", method = "alg2", finder = 1)
{
    rand <- switch(sdist,
                   normal = 1,
                   unif = 2,
                   stop("Selected sampling distribution is not supported!"))
    if(is.null(windows)) {
        windows <- 2**(p+2)
    }
    stopifnot(windows %% 2 == 0)

    pcaoneObj <- switch(method,
                        alg1 = PCAoneAlg1(mat = A, k = k, p = p, q = q, rand = rand, finder),
                        alg2 = PCAoneAlg2(mat = A, k = k, p = p, q = q, rand = rand, windows = windows),
                        stop("Method is not supported!"))
    pcaoneObj$d <- as.vector(pcaoneObj$d)
    pcaoneObj$u <- as.matrix(pcaoneObj$u)
    pcaoneObj$v <- as.matrix(pcaoneObj$v)

    class(pcaoneObj) <- "pcaone"
    return(pcaoneObj)
}
