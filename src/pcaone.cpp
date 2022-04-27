#include <RcppEigen.h>
#if !defined(EIGEN_USE_MKL) // don't use R Lapack.h if MKL is enabled
#include <R_ext/Lapack.h>
#endif

using Rcpp::List;

//' one pass randomized svd in PCAone, adopted from Yu's paer.
//'
//' @param mat the input matrix.
//' @param k top k singular values, k << any(dim(mat)).
//' @param l oversampling.
//' @param p number of power iterations.
//' @export
// [[Rcpp::export]]
List PCAoneH(Eigen::Map<Eigen::MatrixXd> mat, int k, int l, int p) {
  Eigen::Index ncol = mat.cols();
  Eigen::Index nrow = mat.rows();
  Eigen::Index size = l + k;
  Eigen::MatrixXd Omg = Eigen::MatrixXd::Random(ncol, size);
  Eigen::MatrixXd H(ncol, size), G(nrow, size), R(size, size), Rt(size, size);
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

  return List::create(Rcpp::Named("d") = S, Rcpp::Named("u") = U,
                      Rcpp::Named("v") = V);
}
