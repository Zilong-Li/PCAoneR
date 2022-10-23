#include "RsvdEigen.hpp"


//' @title The algorithm1 in PCAone paper.
//' @param mat the input matrix.
//' @param k top k singular values, k << any(dim(mat)).
//' @param p number of power iterations.
//' @param q oversampling.
//' @param rand distribution of random matrix. 1: standard noraml distribution. 2: uniform distribution
//' @param finder method to find othogonal matrix. 1: QR. 2: LU
//' @export
// [[Rcpp::export]]
Rcpp::List PCAoneAlg1(const Eigen::Map<Eigen::MatrixXd> &mat, int k, int p, int q, int rand, int finder = 1) {

  PCAone::RsvdOne<Eigen::MatrixXd> rsvd(mat, k, q, rand);
  rsvd.setRangeFinder(finder);
  rsvd.compute(p);

  return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                            Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                            Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
}

//' @title The algorithm1 in PCAone paper.
//' @param mat the input matrix.
//' @param k top k singular values, k << any(dim(mat)).
//' @param p number of power iterations.
//' @param q oversampling.
//' @param rand distribution of random matrix. 1: standard noraml distribution. 2: uniform distribution
//' @param windows the number of windows. must be a power of 2
//' @export
// [[Rcpp::export]]
Rcpp::List PCAoneAlg2(const Eigen::Map<Eigen::MatrixXd> &mat, int k, int p, int q, int rand, int windows = 64) {

  if (windows % 2 != 0)
    Rcpp::Rcout << "windows %% 2 == 0 has to be met\n";
  if (std::pow(2, p) < windows)
    Rcpp::Rcout << "2^p >= windows has to be met. suggesting p > 6 for windows=64\n";

  PCAone::RsvdOne<Eigen::MatrixXd> rsvd(mat, k, q, rand);
  rsvd.compute(p, windows);

  return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                            Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                            Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
}
