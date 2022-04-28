#include <RcppEigen.h>
#include <random>

using Rcpp::List;

using namespace Eigen;

MatrixXd UniformDist(Eigen::Index n, Eigen::Index m) {
  auto rng = std::default_random_engine {};
  std::uniform_real_distribution<double> dist{0,1};
  const auto uniform{[&](double) { return dist( rng ); }};
  return MatrixXd::NullaryExpr(n, m, uniform);
}

MatrixXd StandardNorm(Eigen::Index n, Eigen::Index m) {
  auto rng = std::default_random_engine {};
  std::normal_distribution<double> dist{0, 1};
  const auto normal{[&](double) { return dist( rng ); }};
  return MatrixXd::NullaryExpr(n, m, normal);
}

//' one pass randomized svd in PCAone, modified from Yu et al 2017.
//' @param mat the input matrix.
//' @param k top k singular values, k << any(dim(mat)).
//' @param l oversampling.
//' @param p number of power iterations.
//' @param rand distribution of random matrix. 1: standard noraml distribution. 2: uniform distribution
//' @export
// [[Rcpp::export]]
List PCAoneY(const Eigen::Map<Eigen::MatrixXd> &mat, int k, int l = 10,
             int p = 3, int rand = 1) {
  Eigen::Index ncol{mat.cols()};
  Eigen::Index nrow{mat.rows()};
  Eigen::Index size = l + k;
  Eigen::MatrixXd Omg(ncol, size), H(ncol, size), G(nrow, size), R(size, size),
      Rt(size, size);
  if (rand == 1) {
    Omg = StandardNorm(ncol, size);
  } else {
    Omg = UniformDist(ncol, size);
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
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(B, Eigen::ComputeThinU |
                                               Eigen::ComputeThinV);
  Eigen::MatrixXd U = G * svd.matrixU().leftCols(k);
  Eigen::MatrixXd V = svd.matrixV().leftCols(k);
  Eigen::VectorXd S = svd.singularValues().head(k);

  return List::create(Rcpp::Named("d") = S,
                      Rcpp::Named("u") = U,
                      Rcpp::Named("v") = V);
}
