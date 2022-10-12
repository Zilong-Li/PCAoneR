#include "RsvdEigen.hpp"


//' @title The one pass randomized svd in PCAone, modified from Yu et al 2017.
//' @param mat the input matrix.
//' @param k top k singular values, k << any(dim(mat)).
//' @param p number of power iterations.
//' @param q oversampling.
//' @param rand distribution of random matrix. 1: standard noraml distribution. 2: uniform distribution
//' @export
// [[Rcpp::export]]
Rcpp::List PCAoneYu(const Eigen::Map<Eigen::MatrixXd> &mat, int k, int p, int q, int rand, int finder) {

  PCAone::RsvdOne<Eigen::MatrixXd> rsvd(mat, k, q, rand);
  rsvd.setRangeFinder(finder);
  rsvd.compute(p);
  // Eigen::VectorXd d = rsvd.singularValues();
  // Rcpp::NumericMatrix u = Rcpp::wrap(rsvd.matrixU());
  // Rcpp::NumericMatrix v = Rcpp::wrap(rsvd.matrixV());

  return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                            Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                            Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
}

//' @title PCAone: the window-based one pass randomized svd, Li et al 2022.
//' @param mat the input matrix.
//' @param k top k singular values, k << any(dim(mat)).
//' @param p number of power iterations.
//' @param q oversampling.
//' @param rand distribution of random matrix. 1: standard noraml distribution. 2: uniform distribution
//' @param windows the number of windows. must be a power of 2
//' @export
// [[Rcpp::export]]
Rcpp::List PCAoneLi(const Eigen::Map<Eigen::MatrixXd> &mat, int k, int p, int q, int rand, int windows) {

  PCAone::RsvdOne<Eigen::MatrixXd> rsvd(mat, k, q, rand);
  rsvd.compute(p, windows);

  return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                            Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                            Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
}
