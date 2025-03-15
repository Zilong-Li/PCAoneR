#include "RsvdEigen.hpp"

// [[Rcpp::depends(RcppEigen)]]

#ifdef EIGEN_USE_MKL_ALL
#include <Eigen/PardisoSupport>
#endif


// [[Rcpp::export]]
Rcpp::List PCAoneAlg1(const Eigen::Map<Eigen::MatrixXd> &mat, int k, int p, int s, int rand) {

  PCAone::RsvdOne<Eigen::MatrixXd> rsvd(mat, k, s, rand);
  int finder = 1;
  rsvd.setRangeFinder(finder);
  rsvd.compute(p);

  return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                            Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                            Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
}

// [[Rcpp::export]]
Rcpp::List PCAoneAlg2(const Eigen::Map<Eigen::MatrixXd> &mat, int k, int p, int s, int rand, int batchs = 64) {

  if (batchs % 2 != 0)
    Rcpp::Rcout << "batchs %% 2 == 0 has to be met\n";
  if (std::pow(2, p) < batchs)
    Rcpp::Rcout << "2^p >= batchs has to be met. suggesting p > 6 for batchs=64\n";

  PCAone::RsvdOne<Eigen::MatrixXd> rsvd(mat, k, s, rand);
  rsvd.compute(p, batchs);

  return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                            Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                            Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
}



// [[Rcpp::export]]
Rcpp::List PCAoneDashSVD(SEXP mat, int k, int p, int s, int rand) {

  Eigen::Map<Eigen::MatrixXd> A = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(mat);
  PCAone::RsvdDash<Eigen::MatrixXd> rsvd(A, k, s, rand);
  rsvd.compute(p);

  return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                            Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                            Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
}

// [[Rcpp::export]]
Rcpp::List DashSVDsparse(SEXP mat, int k, int p, int s, int rand) {

  Eigen::Map<SpMat> A = Rcpp::as<Eigen::Map<SpMat>>(mat);
  PCAone::RsvdDash<SpMat> rsvd(A, k, s, rand);
  rsvd.compute(p);

  return Rcpp::List::create(Rcpp::Named("d") = Rcpp::wrap(rsvd.singularValues()),
                            Rcpp::Named("u") = Rcpp::wrap(rsvd.matrixU()),
                            Rcpp::Named("v") = Rcpp::wrap(rsvd.matrixV()));
}
