#include "RsvdEigen.hpp"


// [[Rcpp::export]]
Rcpp::List PCAoneAlg1(const Eigen::Map<Eigen::MatrixXd> &mat, int k, int p, int q, int rand) {

  PCAone::RsvdOne<Eigen::MatrixXd> rsvd(mat, k, q, rand);
  int finder = 1;
  rsvd.setRangeFinder(finder);
  rsvd.compute(p);

  return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                            Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                            Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
}

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
