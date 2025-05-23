% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcaone.R
\name{pcaone}
\alias{pcaone}
\alias{pcaone.matrix}
\alias{pcaone.dgCMatrix}
\alias{pcaone.dgRMatrix}
\title{Randomized Singular Value Decomposition Algorithm with Window-based Power Iterations from PCAone (Li et al 2022).}
\usage{
pcaone(
  A,
  k = NULL,
  p = 10,
  s = 20,
  method = "auto",
  B = 64,
  shuffle = TRUE,
  opts = list()
)

\method{pcaone}{matrix}(
  A,
  k = NULL,
  p = 10,
  s = 20,
  method = "auto",
  B = 64,
  shuffle = TRUE,
  opts = list()
)

\method{pcaone}{dgCMatrix}(
  A,
  k = NULL,
  p = 10,
  s = 20,
  method = "auto",
  B = 64,
  shuffle = TRUE,
  opts = list()
)

\method{pcaone}{dgRMatrix}(
  A,
  k = NULL,
  p = 7,
  s = 10,
  method = "winsvd",
  B = 64,
  shuffle = TRUE,
  opts = list()
)
}
\arguments{
\item{A}{array_like; \cr
a real/complex \eqn{(m, n)} input matrix (or data frame) to be decomposed.}

\item{k}{integer; \cr
the target rank for the low-rank decomposition. \eqn{k} should satisfy \eqn{k << min(m,n)}.}

\item{p}{integer, optional; \cr
number of power iterations (by default \eqn{p=8}).}

\item{s}{integer, optional; \cr
oversampling parameter (by default \eqn{s=10}).}

\item{method}{string \eqn{c( 'winsvd', 'dashsvd')}, optional; \cr
specifies the different variation of the randomized singular value decomposition : \cr
\eqn{'auto'} (default): automatically choose the method based on size of input matrix. \cr
\eqn{'winsvd'} : window based RSVD in PCAone paper. \cr
\eqn{'dashsvd'} : dynamic shifts based RSVD with  by Feng et al. 2024. \cr}

\item{B}{integer, optional; \cr
the number of batchs for 'winsvd' method. must be a power of 2 (by default \eqn{B=64}).}

\item{shuffle}{logical, optional; \cr
if shuffle the rows of input tall matrix or not for winsvd (by default \eqn{shuffle=TRUE}).}

\item{opts}{list, optional; \cr
options related to PCA, e.g. center and scale. See details}
}
\value{
\code{pcaone} returns a list containing the following three components:
\describe{
\item{d}{  array_like; \cr
singular values; vector of length \eqn{(k)}.
}

\item{u}{  array_like; \cr
left singular vectors; \eqn{(m, k)} or \eqn{(m, nu)} dimensional array.
}

\item{v}{  array_like; \cr
right singular vectors; \eqn{(n, k)} or \eqn{(n, nv)} dimensional array. \cr
}
}
}
\description{
The Randomized Singular Value Decomposition (RSVD) computes the near-optimal low-rank approximation of a rectangular matrix
using a fast probablistic algorithm.
}
\details{
The singular value decomposition (SVD) plays an important role in data analysis, and scientific computing.
Given a rectangular \eqn{(m,n)} matrix \eqn{A}, and a target rank \eqn{k << min(m,n)},
the SVD factors the input matrix \eqn{A} as

\deqn{ A  =  U_{k} diag(d_{k}) V_{k}^\top }{ A = U diag(d) t(V)}

The \eqn{k} left singular vectors are the columns of the
real or complex unitary matrix \eqn{U}. The \eqn{k} right singular vectors are the columns
of the real or complex unitary matrix \eqn{V}. The \eqn{k} dominant singular values are the
entries of \eqn{d}, and non-negative and real numbers.

\eqn{q} is an oversampling parameter to improve the approximation.
A value of at least 10 is recommended, and \eqn{q=10} is set by default.

The parameter \eqn{p} specifies the number of power (subspace) iterations
to reduce the approximation error. The power scheme is recommended,
especially when the singular values decay slowly. However, computing power iterations increases the
computational costs. Even though most RSVD implementations recommend \eqn{p=3} power iterations by default,
it's always sufficient to run only few power iterations where our window-based power iterations (\eqn{'alg2'})
come to play. We recommend using \eqn{windows=64} and \eqn{p>=7} for pcaone algorithm2. As it is designed for large dataset,
we recommend using \eqn{'alg2'} when \eqn{max(n,m) > 5000}.

If \eqn{k > (min(n,m)/4)}, a deterministic partial or truncated \code{\link{svd}}
algorithm might be faster.
}
\note{
The singular vectors are not unique and only defined up to sign.
If a left singular vector has its sign changed, changing the sign of the corresponding right vector
gives an equivalent decomposition.
}
\examples{
library('pcaone')
load(system.file("extdata", "popgen.rda", package="pcaone") )
A <- popgen - rowMeans(popgen)
res <- pcaone(A, k = 40)
str(res)
}
\references{
\itemize{
\item Z. Li, J Meisner, A Albrechtsen. "Fast and accurate out-of-core PCA framework for large scale biobank data" (2023)
\doi{10.1101/gr.277525.122}.
}
}
\author{
Zilong Li \email{zilong.dk@gmail.com}
}
