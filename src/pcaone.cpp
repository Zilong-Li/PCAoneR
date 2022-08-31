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
List PCAoneY(const Eigen::Map<Eigen::MatrixXd> &mat, int k, int p = 3,
             int q = 10, int rand = 1) {
  Eigen::Index ncol{mat.cols()};
  Eigen::Index nrow{mat.rows()};
  Eigen::Index size = q + k;
  Eigen::MatrixXd Omg(ncol, size), H(ncol, size), G(nrow, size), R(size, size),
      Rt(size, size);
  auto rng = std::default_random_engine{};
  if (rand == 1) {
    Omg = PCAone::StandardNormalRandom<Eigen::MatrixXd, std::default_random_engine>(ncol, size, rng);
  } else {
    Omg = PCAone::UniformRandom<Eigen::MatrixXd, std::default_random_engine>(ncol, size, rng);
  }
  G.noalias() = mat * Omg;
  H.noalias() = mat.transpose() * G;
  if (p > 0) {
    for (int i = 0; i < p; ++i) {
      Eigen::HouseholderQR<Eigen::Ref<Eigen::MatrixXd>> qr(H);
      H.noalias() = qr.householderQ() * Eigen::MatrixXd::Identity(ncol, size);
      G.noalias() = mat * H;
      H.noalias() = mat.transpose() * G;
    }
  }
  {
    Eigen::HouseholderQR<Eigen::Ref<Eigen::MatrixXd>> qr(G);
    R.noalias() = Eigen::MatrixXd::Identity(size, nrow) *
                  qr.matrixQR().template triangularView<Eigen::Upper>();
    G.noalias() = qr.householderQ() * Eigen::MatrixXd::Identity(nrow, size);
  }
  {
    Eigen::HouseholderQR<Eigen::Ref<Eigen::MatrixXd>> qr(G);
    Rt.noalias() = Eigen::MatrixXd::Identity(size, nrow) *
                   qr.matrixQR().template triangularView<Eigen::Upper>();
    G.noalias() = qr.householderQ() * Eigen::MatrixXd::Identity(nrow, size);
  }

  R = Rt * R;
  Eigen::MatrixXd B = R.transpose().colPivHouseholderQr().solve(H.transpose());
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(B, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::MatrixXd U = G * svd.matrixU().leftCols(k);
  Eigen::MatrixXd V = svd.matrixV().leftCols(k);
  Eigen::VectorXd S = svd.singularValues().head(k);

  return List::create(Rcpp::Named("d") = S,
                      Rcpp::Named("u") = U,
                      Rcpp::Named("v") = V);
}
