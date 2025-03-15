#include "RsvdEigen.hpp"

// [[Rcpp::depends(RcppEigen)]]

#ifdef EIGEN_USE_MKL_ALL
#include <Eigen/PardisoSupport>
#endif

using SpMat = Eigen::SparseMatrix<double, Eigen::ColMajor>;
using SpRMat = Eigen::SparseMatrix<double, Eigen::RowMajor>;


// [[Rcpp::export]]
Rcpp::List ssvd(SEXP mat, int k, int p, int s, int rand) {

  Eigen::Map<Eigen::MatrixXd> A = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(mat);
  PCAone::RsvdOne<Eigen::MatrixXd> rsvd(A, k, s, rand);
  int finder = 1;
  rsvd.setRangeFinder(finder);
  rsvd.compute(p);

  return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                            Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                            Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
}

// [[Rcpp::export]]
Rcpp::List winsvd(SEXP mat, int k, int p, int s, int rand, int batchs = 64) {

  if (batchs % 2 != 0)
    Rcpp::Rcout << "batchs %% 2 == 0 has to be met\n";
  if (std::pow(2, p) < batchs)
    Rcpp::Rcout << "2^p >= batchs has to be met. suggesting p > 6 for batchs=64\n";

  Eigen::Map<Eigen::MatrixXd> A = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(mat);
  PCAone::RsvdOne<Eigen::MatrixXd> rsvd(A, k, s, rand);
  rsvd.compute(p, batchs);

  return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                            Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                            Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
}



// [[Rcpp::export]]
Rcpp::List dashsvd(SEXP mat, int k, int p, int s, int rand) {

  Eigen::Map<Eigen::MatrixXd> A = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(mat);
  PCAone::RsvdDash<Eigen::MatrixXd> rsvd(A, k, s, rand);
  rsvd.compute(p);

  return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                            Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                            Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
}

// [[Rcpp::export]]
Rcpp::List dashsvd_sparse_col(SEXP mat, int k, int p, int s, int rand) {

  Eigen::Map<SpMat> A = Rcpp::as<Eigen::Map<SpMat>>(mat);
  PCAone::RsvdDash<SpMat> rsvd(A, k, s, rand);
  rsvd.compute(p);

  return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                            Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                            Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
}

// [[Rcpp::export]]
Rcpp::List dashsvd_sparse_row(SEXP mat, int k, int p, int s, int rand) {

  Eigen::Map<SpRMat> A = Rcpp::as<Eigen::Map<SpRMat>>(mat);
  PCAone::RsvdDash<SpRMat> rsvd(A, k, s, rand);
  rsvd.compute(p);

  return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                            Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                            Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
}
