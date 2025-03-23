/***********************************************************************
 * @file        RsvdEigen.hpp
 * @author      Zilong Li
 * @email       zilong.dk@gmail.com
 * @brief       Window-based Randomized Singular Value Decomposition
 * Copyright (C) 2021-2025 Zilong Li
 ***********************************************************************/

#ifndef RSVDEIGEN_H_
#define RSVDEIGEN_H_

#include "Rand.hpp"
#include <RcppEigen.h>
#include <memory>
#include <stdexcept>

namespace PCAone
{

template<typename MatrixType>
void flipOmg(MatrixType & Omg2, MatrixType & Omg)
{
  for(Eigen::Index i = 0; i < Omg.cols(); ++i)
  {
    // if signs of half of values are flipped then correct signs.
    if((Omg2.col(i) - Omg.col(i)).array().abs().sum() > 2 * (Omg2.col(i) + Omg.col(i)).array().abs().sum())
    {
      Omg.col(i) *= -1;
    }
  }
  Omg2 = Omg;
}

template<typename MatrixType>
class RsvdOne
{
private:
  using Index = Eigen::Index;
  using Scalar = typename MatrixType::Scalar;
  using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>; // dense matrix
  using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
  using ConstRefVec = const Eigen::Ref<const Vector>;
  using ConstGenericMatrix = const Eigen::Ref<const MatrixType>;

  ConstGenericMatrix A;
  const int k, os, rand, size;
  const bool byrow;
  ConstRefVec center;
  ConstRefVec scale;
  uint32_t batch_size;
  
  bool dopca = false; // if set center and scale, then do pca
  bool trans; //  work on A.T
  Index nrow, ncol;  // can change
  Matrix Omg; // Omega to update
  Matrix Ab;  // block of normalized matrix we work on for PCA
  Matrix b_leftSingularVectors;
  Vector b_singularValues;
  Matrix b_rightSingularVectors;

public:
  RsvdOne(ConstGenericMatrix & A_, int k_, int os_, int rand_, bool byrow_, ConstRefVec & ctr, ConstRefVec& scl)
  : A(A_), k(k_), os(os_), rand(rand_), size(k_+os_), byrow(byrow_), center(ctr), scale(scl)
  {
    trans = A.rows() >= A.cols() ? false : true;
    if(center.size() > 1 && scale.size() > 1) dopca = true;
    if(dopca) trans = byrow ? false : true;

    // whether we should work on the transposed matrix
    if(trans == false)
    {
      nrow = A.rows();
      ncol = A.cols();
    }
    else
    {
      nrow = A.cols();
      ncol = A.rows();
    }

    auto randomEngine = std::default_random_engine{};
    if(rand == 1)
      Omg = StandardNormalRandom<Matrix, std::default_random_engine>(ncol, size, randomEngine);
    else
      Omg = UniformRandom<Matrix, std::default_random_engine>(ncol, size, randomEngine);
  }

  ~RsvdOne() {}

  void compute(int p, int batches)
  {
    Matrix H(ncol, size), G(nrow, size), R(size, size), Rt(size, size);
    // G = A * Omega; H = A.transpose() * G;
    computeGandH(G, H, p, batches);

    {
      Eigen::HouseholderQR<Eigen::Ref<Matrix>> qr(G);
      R.noalias() = Matrix::Identity(size, nrow) * qr.matrixQR().template triangularView<Eigen::Upper>();
      G.noalias() = qr.householderQ() * Matrix::Identity(nrow, size);
    }
    {
      Eigen::HouseholderQR<Eigen::Ref<Matrix>> qr(G);
      Rt.noalias() = Matrix::Identity(size, nrow) * qr.matrixQR().template triangularView<Eigen::Upper>();
      G.noalias() = qr.householderQ() * Matrix::Identity(nrow, size);
    }

    R = Rt * R;

    // R.T * B = H.T => lapack dtrtrs()
    Matrix B = R.transpose().fullPivHouseholderQr().solve(H.transpose());
    Eigen::JacobiSVD<Matrix> svd(B, Eigen::ComputeThinU | Eigen::ComputeThinV);
    b_leftSingularVectors.noalias() = G * svd.matrixU().leftCols(k);
    b_rightSingularVectors = svd.matrixV().leftCols(k);
    b_singularValues = svd.singularValues().head(k);
  }
  
  inline Matrix matrixU() const
  {
    if(trans)
      return b_rightSingularVectors;
    else
      return b_leftSingularVectors;
  }

  inline Matrix matrixV() const
  {
    if(trans)
      return b_leftSingularVectors;
    else
      return b_rightSingularVectors;
  }

  inline Vector singularValues() const
  {
    return b_singularValues;
  }

private:

  void computeGandH(Matrix & G, Matrix & H, int p, int batches)
  {
    if(batches % 2 != 0) throw std::runtime_error("batches must be a power of 2, i.e. batches=2^x.\n");
    if(std::pow(2, p) < batches) throw std::runtime_error("pow(2, p) >= batches is not true\n");
    // nrow of a tall matrix is used for sliding window 
    batch_size = (uint32_t)std::ceil((double)nrow / batches);
    if(batch_size < batches) throw std::runtime_error("window size is smaller than number of batches \n");

    uint32_t start, stop, nb;
    Matrix H1 = Matrix::Zero(ncol, size);
    Matrix H2 = Matrix::Zero(ncol, size);
    Matrix Omg2;
    int band = 1;
    for(uint32_t pi = 0; pi <= p; pi++)
    {
      if(pi == 0) Omg2 = Omg;
      if(std::pow(2, pi) >= batches)
      {
        // reset H1, H2 to zero
        H1.setZero();
        H2.setZero();
      }
      band = std::fmin(band * 2, batches);
      for(int b = 0, i = 1, j = 1; b < batches; ++b, ++i, ++j)
      {
        start = b * batch_size;
        stop = (b + 1) * batch_size >= nrow ? nrow - 1 : (b + 1) * batch_size - 1;
        nb = stop - start + 1;
        if(!dopca)
        {
          if(trans)
          {
            G.middleRows(start, nb).noalias() = A.middleCols(start, nb).transpose() * Omg;
            if(i <= band / 2)
              H1.noalias() += A.middleCols(start, nb) * G.middleRows(start, nb);
            else
              H2.noalias() += A.middleCols(start, nb) * G.middleRows(start, nb);
          }
          else
          {
            G.middleRows(start, nb).noalias() = A.middleRows(start, nb) * Omg;
            if(i <= band / 2)
              H1.noalias() += A.middleRows(start, nb).transpose() * G.middleRows(start, nb);
            else
              H2.noalias() += A.middleRows(start, nb).transpose() * G.middleRows(start, nb);
          }
        }
        else
        {
          // here we can only trans it if byrow is false
          if(trans)
          {
            center_and_scale(Ab, A.middleCols(start, nb), center, scale, false, start, nb);
            G.middleRows(start, nb).noalias() = Ab.transpose() * Omg;
            if(i <= band / 2)
              H1.noalias() += Ab * G.middleRows(start, nb);
            else
              H2.noalias() += Ab * G.middleRows(start, nb);
          }
          else
          {
            center_and_scale(Ab, A.middleRows(start, nb), center, scale, true, start, nb);
            G.middleRows(start, nb).noalias() = Ab * Omg;
            if(i <= band / 2)
              H1.noalias() += Ab.transpose() * G.middleRows(start, nb);
            else
              H2.noalias() += Ab.transpose() * G.middleRows(start, nb);
          }
        }

        // use the first quarter band of succesive iteration (H1)
        // for extra power iteration updates with the last used band (H2)
        const bool adjacent = (pi > 0 && (b + 1) == std::pow(2, pi - 1) && std::pow(2, pi) < batches);
        if((b + 1) < band && !adjacent) continue;
        if(!((i == band) || (i == band / 2) || adjacent)) continue;
        H = H1 + H2;
        Eigen::HouseholderQR<Matrix> qr(H);
        Omg.noalias() = qr.householderQ() * Matrix::Identity(ncol, size);
        flipOmg(Omg2, Omg);
        if(i == band)
        {
          H1.setZero();
          i = 0;
        }
        else
        {
          H2.setZero();
        }
      }
    }
  }


};

template<typename MatrixType>
class RsvdDash
{
private:
  using Index = Eigen::Index;
  using Scalar = typename MatrixType::Scalar;
  using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
  using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
  using MapConstVec = const Eigen::Ref<const Vector>;
  using ConstGenericMatrix = const Eigen::Ref<const MatrixType>;
  
  ConstGenericMatrix A;
  const int k, os, size, rand;
  const Index nrow, ncol;
  bool dopca = false;
  const bool byrow;
  bool tall;
  MapConstVec center;
  MapConstVec scale;
  Matrix Omg, Ab;
  Matrix b_leftSingularVectors;
  Vector b_singularValues;
  Matrix b_rightSingularVectors;
  
public:
  RsvdDash(ConstGenericMatrix & A_, int k_, int os_, int rand_, bool byrow_, MapConstVec & ctr, MapConstVec& scl)
  : A(A_), k(k_), os(os_), size(k_ + os_), rand(rand_), nrow(A_.rows()), ncol(A_.cols()), byrow(byrow_), center(ctr), scale(scl)
  {
    if(size > nrow || size > ncol)
      throw std::runtime_error("k(pc) + s(oversampling) must be not greater than the dimension of data");
    tall = nrow > ncol ? true : false;
    auto randomEngine = std::default_random_engine{};
    if(rand == 1)
    {
      Omg = StandardNormalRandom<Matrix, std::default_random_engine>(tall ? nrow : ncol, size, randomEngine);
    }
    else
    {
      Omg = UniformRandom<Matrix, std::default_random_engine>(tall ? nrow : ncol, size, randomEngine);
    }
    b_leftSingularVectors = Matrix::Zero(nrow, size);
    b_rightSingularVectors = Matrix::Zero(ncol, size);
    b_singularValues = Vector::Zero(size);
    if(center.size() > 1 && scale.size() > 1) {
      dopca = true;
      Ab = center_and_scale(A, center, scale, byrow);
    }
  }

  ~RsvdDash() {}

  inline void compute(uint32_t p)
  {
    if(tall)
    {
      compute_tall(p);
    }
    else
    {
      compute_wide(p);
    }
  }

  void compute_tall(uint32_t p)
  {
    if(dopca) {
      Omg = Ab.transpose() * Omg;
    } else {
      Omg = A.transpose() * Omg;
    }
    eigSVD(Omg, b_rightSingularVectors, b_singularValues, b_leftSingularVectors);
    double alpha = 0.0;
    for(uint32_t i = 0; i < p; i++)
    {
      if(dopca) {
        Omg.noalias() = Ab.transpose() * (Ab * b_rightSingularVectors) - alpha * b_rightSingularVectors;
      } else {
        Omg.noalias() = A.transpose() * (A * b_rightSingularVectors) - alpha * b_rightSingularVectors;
      }
      eigSVD(Omg, b_rightSingularVectors, b_singularValues, b_leftSingularVectors);
      if(b_singularValues[size - 1] > alpha) alpha = (alpha + b_singularValues[size - 1]) / 2.0;
    }
    if(dopca) {
      eigSVD(Ab * b_rightSingularVectors, b_leftSingularVectors, b_singularValues, Omg);
    } else {
      eigSVD(A * b_rightSingularVectors, b_leftSingularVectors, b_singularValues, Omg);
    }
    b_rightSingularVectors = b_rightSingularVectors * Omg;
  }

  void compute_wide(uint32_t p)
  {
    if(dopca) {
      Omg = Ab * Omg;
    } else {
      Omg = A * Omg;
    }
    eigSVD(Omg, b_leftSingularVectors, b_singularValues, b_rightSingularVectors);
    double alpha = 0.0;
    for(uint32_t i = 0; i < p; i++)
    {
      if(dopca) {
        Omg.noalias() = Ab * (Ab.transpose() * b_leftSingularVectors) - alpha * b_leftSingularVectors;
      } else {
        Omg.noalias() = A * (A.transpose() * b_leftSingularVectors) - alpha * b_leftSingularVectors;
      }
      eigSVD(Omg, b_leftSingularVectors, b_singularValues, b_rightSingularVectors);
      if(b_singularValues[size - 1] > alpha) alpha = (alpha + b_singularValues[size - 1]) / 2.0;
    }
    if(dopca) {
      eigSVD(Ab.transpose() * b_leftSingularVectors, b_rightSingularVectors, b_singularValues, Omg);
    } else {
      eigSVD(A.transpose() * b_leftSingularVectors, b_rightSingularVectors, b_singularValues, Omg);
    }
    b_leftSingularVectors = b_leftSingularVectors * Omg;
  }

  /// A is tall and thin
  void eigSVD(const Matrix & mat, Matrix & U, Vector & S, Matrix & V)
  {
    Matrix C = mat.transpose() * mat;
    Eigen::SelfAdjointEigenSolver<Matrix> eigen_solver(C);
    if(eigen_solver.info() != Eigen::Success)
    {
      throw std::runtime_error("Eigendecomposition failed");
    }
    // Reverse order for descending eigenvalues
    V.noalias() = eigen_solver.eigenvectors().rowwise().reverse();
    Vector eigenvalues = eigen_solver.eigenvalues().reverse();

    // Compute singular values (sqrt of eigenvalues)
    S.noalias() = eigenvalues.cwiseSqrt();

    // Handle numerical stability for singular values
    const double tolerance = 1e-10 * S.maxCoeff() * mat.cols();
    Vector inv_S = (S.array().abs() > tolerance).select(S.array().inverse(), 0.0);

    U.noalias() = mat * V * inv_S.asDiagonal();
  }

  // return a copy or maybe a reference ?
  inline Matrix matrixU() const
  {
    return b_leftSingularVectors.leftCols(k);
  }

  inline Matrix matrixV() const
  {
    return b_rightSingularVectors.leftCols(k);
  }

  inline Vector singularValues() const
  {
    return b_singularValues.head(k);
  }
};

} // namespace PCAone

#endif // RSVDEIGEN_H_
