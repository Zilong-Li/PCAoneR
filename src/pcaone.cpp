#include "RsvdEigen.hpp"

using Rcpp::List;

using namespace Eigen;

//' @title The one pass randomized svd in PCAone, modified from Yu et al 2017.
//' @param mat the input matrix.
//' @param k top k singular values, k << any(dim(mat)).
//' @param p number of power iterations.
//' @param q oversampling.
//' @param rand distribution of random matrix. 1: standard noraml distribution. 2: uniform distribution
//' @export
// [[Rcpp::export]]
List PCAoneYu(const Eigen::Map<Eigen::MatrixXd> &mat, int k, int p, int q, int rand) {

  PCAone::RsvdOne<Eigen::MatrixXd> rsvd(mat, k, q, rand);
  rsvd.compute(p);

  return List::create(Rcpp::Named("d") = rsvd.singularValues(),
                      Rcpp::Named("u") = rsvd.matrixU(),
                      Rcpp::Named("v") = rsvd.matrixV());
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
List PCAoneLi(const Eigen::Map<Eigen::MatrixXd> &mat, int k, int p, int q, int rand, int windows) {

  PCAone::RsvdOne<Eigen::MatrixXd> rsvd(mat, k, q, rand);
  rsvd.compute(p, windows);

  return List::create(Rcpp::Named("d") = rsvd.singularValues(),
                      Rcpp::Named("u") = rsvd.matrixU(),
                      Rcpp::Named("v") = rsvd.matrixV());
}
