#include "RsvdEigen.hpp"

// [[Rcpp::depends(RcppEigen)]]

#ifdef EIGEN_USE_MKL_ALL
#include <Eigen/PardisoSupport>
#endif

using SpMat = Eigen::SparseMatrix<double, Eigen::ColMajor>;
using SpRMat = Eigen::SparseMatrix<double, Eigen::RowMajor>;

// [[Rcpp::export]]
Rcpp::List svd_dense(SEXP mat, int k, int p, int s, int batchs, Rcpp::List params_pca) {
  Eigen::Map<Eigen::MatrixXd> A = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(mat);
  Eigen::Map<Eigen::VectorXd> center = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(params_pca["center"]);
  Eigen::Map<Eigen::VectorXd> scale = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(params_pca["scale"]);

  int  rand  = Rcpp::as<int>(params_pca["rand"]);
  int  method= Rcpp::as<int>(params_pca["method"]);
  bool dopca = Rcpp::as<bool>(params_pca["dopca"]);
  bool byrow = Rcpp::as<bool>(params_pca["byrow"]);

  if(method == 3) {
    PCAone::RsvdDash<Eigen::MatrixXd> rsvd(A, k, s, rand, byrow, center, scale);
    rsvd.compute(p);
    return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                              Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                              Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
  }
  
  if(method == 2) batchs = 0;  // ssvd
  PCAone::RsvdOne<Eigen::MatrixXd> rsvd(A, k, s, rand, byrow, center, scale);
  rsvd.compute(p, batchs);

  return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                            Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                            Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
}


// [[Rcpp::export]]
Rcpp::List svd_sparse_col(SEXP mat, int k, int p, int s, int batchs, Rcpp::List params_pca) {
  Eigen::Map<SpMat> A = Rcpp::as<Eigen::Map<SpMat>>(mat);
  Eigen::Map<Eigen::VectorXd> center = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(params_pca["center"]);
  Eigen::Map<Eigen::VectorXd> scale = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(params_pca["scale"]);

  int  rand  = Rcpp::as<int>(params_pca["rand"]);
  int  method= Rcpp::as<int>(params_pca["method"]);
  bool dopca = Rcpp::as<bool>(params_pca["dopca"]);
  bool byrow = Rcpp::as<bool>(params_pca["byrow"]);

  if(method == 3) {
    PCAone::RsvdDash<SpMat> rsvd(A, k, s, rand, byrow, center, scale);
    rsvd.compute(p);
    return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                              Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                              Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
  }
  
  if(method == 2) batchs = 0;  // ssvd
  PCAone::RsvdOne<SpMat> rsvd(A, k, s, rand, byrow, center, scale);
  rsvd.compute(p, batchs);

  return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                            Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                            Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
}

// [[Rcpp::export]]
Rcpp::List svd_sparse_row(Rcpp::S4 mat, int k, int p, int s, int batchs, Rcpp::List params_pca) {
  // Extract dimensions (nrow, ncol)
  Rcpp::IntegerVector dims = mat.slot("Dim");
  int nrow = dims[0];
  int ncol = dims[1];
  // Extract slots
  Rcpp::IntegerVector q = mat.slot("p"); // Row pointers 
  Rcpp::IntegerVector j = mat.slot("j"); // Column indices (0-based)
  Rcpp::NumericVector x = mat.slot("x"); // Non-zero values
  // Manually map to Eigen::SparseMatrix<double, RowMajor>
  Eigen::Map<SpRMat> A(nrow, ncol,         // Explicitly set dimensions 
                       x.size(),           // Number of non-zero elements
                       q.begin(),          // Row pointers (size = nrow + 1)
                       j.begin(),          // Column indices (0-based)
                       x.begin()           // Values
                       );
  // Extract pca opts
  int  rand  = Rcpp::as<int>(params_pca["rand"]);
  int  method= Rcpp::as<int>(params_pca["method"]);
  bool dopca = Rcpp::as<bool>(params_pca["dopca"]);
  bool byrow = Rcpp::as<bool>(params_pca["byrow"]);
  Eigen::Map<Eigen::VectorXd> center = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(params_pca["center"]);
  Eigen::Map<Eigen::VectorXd> scale = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(params_pca["scale"]);

  if(method == 3) {
    PCAone::RsvdDash<SpRMat> rsvd(A, k, s, rand, byrow, center, scale);
    rsvd.compute(p);
    return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                              Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                              Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
  }
  
  if(method == 2) batchs = 0;  // ssvd
  PCAone::RsvdOne<SpRMat> rsvd(A, k, s, rand, byrow, center, scale);
  rsvd.compute(p, batchs);

  return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                            Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                            Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
}


