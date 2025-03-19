#include "RsvdEigen.hpp"

// [[Rcpp::depends(RcppEigen)]]

#ifdef EIGEN_USE_MKL_ALL
#include <Eigen/PardisoSupport>
#endif

using SpMat = Eigen::SparseMatrix<double, Eigen::ColMajor>;
using SpRMat = Eigen::SparseMatrix<double, Eigen::RowMajor>;

// [[Rcpp::export]]
Rcpp::List winsvd(SEXP mat, int k, int p, int s, int rand, int batchs, Rcpp::List params_pca) {
  bool dopca = Rcpp::as<bool>(params_pca["dopca"]);
  bool byrow = Rcpp::as<bool>(params_pca["byrow"]);
  Rcpp::NumericVector center  = params_pca["center"];
  Rcpp::NumericVector scale  = params_pca["scale"];

  Eigen::Map<Eigen::MatrixXd> A = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(mat);
  PCAone::RsvdOne<Eigen::MatrixXd> rsvd(A, k, s, rand);
  if(dopca) rsvd.setCenterScale(center.length(), center.begin(), scale.begin(), byrow);
  rsvd.compute(p, batchs);

  return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                            Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                            Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
}

// [[Rcpp::export]]
Rcpp::List winsvd_sparse_col(SEXP mat, int k, int p, int s, int rand, int batchs) {
  Eigen::Map<SpMat> A = Rcpp::as<Eigen::Map<SpMat>>(mat);
  PCAone::RsvdOne<SpMat> rsvd(A, k, s, rand);
  rsvd.compute(p, batchs);

  return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                            Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                            Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
}

// [[Rcpp::export]]
Rcpp::List winsvd_sparse_row(Rcpp::S4 mat, int k, int p, int s, int rand, int batchs) {
    // Extract dimensions (nrow, ncol)
  Rcpp::IntegerVector dims = mat.slot("Dim");
  int nrow = dims[0];
  int ncol = dims[1];

  // Extract slots
  Rcpp::IntegerVector q = mat.slot("p"); // Row pointers 
  Rcpp::IntegerVector j = mat.slot("j"); // Column indices (0-based)
  Rcpp::NumericVector x = mat.slot("x"); // Non-zero values
  // Manually map to Eigen::SparseMatrix<double, RowMajor>
  Eigen::Map<SpRMat> A(
    nrow, ncol,         // Explicitly set dimensions 
    x.size(),           // Number of non-zero elements
    q.begin(),          // Row pointers (size = nrow + 1)
    j.begin(),          // Column indices (0-based)
    x.begin()           // Values
  );

  PCAone::RsvdOne<SpRMat> rsvd(A, k, s, rand);
  rsvd.compute(p, batchs);

  return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                            Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                            Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
}

// [[Rcpp::export]]
Rcpp::List ssvd(SEXP mat, int k, int p, int s, int rand, Rcpp::List params_pca) {
  bool dopca = Rcpp::as<bool>(params_pca["dopca"]);
  bool byrow = Rcpp::as<bool>(params_pca["byrow"]);
  Rcpp::NumericVector center  = params_pca["center"];
  Rcpp::NumericVector scale  = params_pca["scale"];
  
  Eigen::Map<Eigen::MatrixXd> A = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(mat);
  PCAone::RsvdOne<Eigen::MatrixXd> rsvd(A, k, s, rand);
  if(dopca) rsvd.setCenterScale(center.length(), center.begin(), scale.begin(), byrow);
  rsvd.compute(p);

  return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                            Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                            Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
}

// [[Rcpp::export]]
Rcpp::List ssvd_sparse_col(SEXP mat, int k, int p, int s, int rand) {
  Eigen::Map<SpMat> A = Rcpp::as<Eigen::Map<SpMat>>(mat);
  PCAone::RsvdOne<SpMat> rsvd(A, k, s, rand);
  rsvd.compute(p);

  return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                            Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                            Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
}

// [[Rcpp::export]]
Rcpp::List ssvd_sparse_row(Rcpp::S4 mat, int k, int p, int s, int rand) {
  // Eigen::Map<SpRMat> A = Rcpp::as<Eigen::Map<SpRMat>>(mat);
  // Extract dimensions (nrow, ncol)
  Rcpp::IntegerVector dims = mat.slot("Dim");
  int nrow = dims[0];
  int ncol = dims[1];

  // Extract slots
  Rcpp::IntegerVector q = mat.slot("p"); // Row pointers 
  Rcpp::IntegerVector j = mat.slot("j"); // Column indices (0-based)
  Rcpp::NumericVector x = mat.slot("x"); // Non-zero values
  // Manually map to Eigen::SparseMatrix<double, RowMajor>
  Eigen::Map<SpRMat> A(
    nrow, ncol,         // Explicitly set dimensions 
    x.size(),           // Number of non-zero elements
    q.begin(),          // Row pointers (size = nrow + 1)
    j.begin(),          // Column indices (0-based)
    x.begin()           // Values
  );

  PCAone::RsvdOne<SpRMat> rsvd(A, k, s, rand);
  rsvd.compute(p);
  return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                            Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                            Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
}

// [[Rcpp::export]]
Rcpp::List dashsvd(SEXP mat, int k, int p, int s, int rand, Rcpp::List params_pca) {

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
Rcpp::List dashsvd_sparse_row(Rcpp::S4 mat, int k, int p, int s, int rand) {
  // Extract dimensions (nrow, ncol)
  Rcpp::IntegerVector dims = mat.slot("Dim");
  int nrow = dims[0];
  int ncol = dims[1];

  // Extract slots
  Rcpp::IntegerVector q = mat.slot("p"); // Row pointers 
  Rcpp::IntegerVector j = mat.slot("j"); // Column indices (0-based)
  Rcpp::NumericVector x = mat.slot("x"); // Non-zero values
  // Manually map to Eigen::SparseMatrix<double, RowMajor>
  Eigen::Map<SpRMat> A(
    nrow, ncol,         // Explicitly set dimensions 
    x.size(),           // Number of non-zero elements
    q.begin(),          // Row pointers (size = nrow + 1)
    j.begin(),          // Column indices (0-based)
    x.begin()           // Values
  );

  PCAone::RsvdDash<SpRMat> rsvd(A, k, s, rand);
  rsvd.compute(p);

  return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                            Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                            Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
}
