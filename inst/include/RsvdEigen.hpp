/***********************************************************************
 * @file        RsvdEigen.hpp
 * @author      Zilong Li
 * @email       zilong.dk@gmail.com
 * @brief       Various Randomized Singular Value Decomposition in Eigen
 * Copyright (C) 2021 Zilong Li
 ***********************************************************************/

#ifndef RSVDEIGEN_H_
#define RSVDEIGEN_H_

#include "Rand.hpp"
#include <stdexcept>

namespace PCAone
{

template<typename MatrixType>
void flipOmg(MatrixType & Omg2, MatrixType & Omg)
{
    for(Eigen::Index i = 0; i < Omg.cols(); ++i)
    {
        // if signs of half of values are flipped then correct signs.
        if((Omg2.col(i) - Omg.col(i)).array().abs().sum()
           > 2 * (Omg2.col(i) + Omg.col(i)).array().abs().sum())
        {
            Omg.col(i) *= -1;
        }
    }
    Omg2 = Omg;
}

// RsvdOp: rows, cols, ranks, oversamples, computeGandH()
template<typename MatrixType>
class RsvdOpOnePass
{
  private:
    using ConstGenericMatrix = const Eigen::Ref<const MatrixType>;
    using Index = Eigen::Index;
    using Scalar = typename MatrixType::Scalar;
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>; // dense matrix
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    ConstGenericMatrix mat;
    Index nrow, ncol;
    Matrix Omg;
    int k, os, size, rand;
    bool trans; // if matrix is wide then flip the matrix dimension

  public:
    int finder = 1;
    RsvdOpOnePass(ConstGenericMatrix & mat_, int k_, int os_ = 10, int rand_ = 1)
    : mat(mat_), k(k_), os(os_), size(k_ + os_), rand(rand_)
    {
        if(mat.rows() >= mat.cols())
        {
            trans = false;
            nrow = mat.rows();
            ncol = mat.cols();
        }
        else
        {
            trans = true;
            nrow = mat.cols();
            ncol = mat.rows();
        }

        auto randomEngine = std::default_random_engine{};
        if(rand == 1)
        {
            Omg = StandardNormalRandom<Matrix, std::default_random_engine>(ncol, size, randomEngine);
        }
        else
        {
            Omg = UniformRandom<Matrix, std::default_random_engine>(ncol, size, randomEngine);
        }
    }

    virtual ~RsvdOpOnePass() {}

    Index rows() const
    {
        return nrow;
    }
    Index cols() const
    {
        return ncol;
    }
    Index ranks() const
    {
        return k;
    }
    Index oversamples() const
    {
        return os;
    }

    void computeGandH(Matrix & G, Matrix & H, uint32_t p)
    {
        if(trans)
        {
            G.noalias() = mat.transpose() * Omg;
            H.noalias() = mat * G;
        }
        else
        {
            G.noalias() = mat * Omg;
            H.noalias() = mat.transpose() * G;
        }
        if(p > 0)
        {
            for(uint32_t pi = 0; pi < p; pi++)
            {
                if(finder == 1)
                {
                    Eigen::HouseholderQR<Eigen::Ref<Matrix>> qr(H);
                    Omg.noalias() = qr.householderQ() * Matrix::Identity(cols(), size);
                }
                else if(finder == 2)
                {
                    Eigen::FullPivLU<Eigen::Ref<Matrix>> lu(H);
                    Omg.setIdentity(cols(), size);
                    Omg.template triangularView<Eigen::StrictlyLower>() = lu.matrixLU();
                }
                else
                {
                    throw std::invalid_argument("finder must be 1 or 2");
                }
                if(trans)
                {
                    G.noalias() = mat.transpose() * Omg;
                    H.noalias() = mat * G;
                }
                else
                {
                    G.noalias() = mat * Omg;
                    H.noalias() = mat.transpose() * G;
                }
            }
        }
    }

    void computeGandH(Matrix & G, Matrix & H, uint32_t p, uint32_t batchs)
    {
        if(batchs % 2 != 0) throw std::runtime_error("batchs must be a power of 2, ie. batchs=2^x.\n");
        if(std::pow(2, p) < batchs) throw std::runtime_error("pow(2, p) >= batchs is not true\n");
        uint32_t blocksize = (unsigned int)std::ceil((double)nrow / batchs);
        if(blocksize < batchs)
            throw std::runtime_error(
                "window size is smaller than number of batchs because given matrix is too small. please "
                "consider other methods or adjust parameter batchs.\n");
        if(trans)
        {
            G.noalias() = mat.transpose() * Omg;
            H.noalias() = mat * G;
        }
        else
        {
            G.noalias() = mat * Omg;
            H.noalias() = mat.transpose() * G;
        }

        uint32_t start_idx, stop_idx, actual_block_size;
        Matrix H1 = Matrix::Zero(ncol, size);
        Matrix H2 = Matrix::Zero(ncol, size);
        Matrix Omg2;
        uint32_t band = 1;
        for(uint32_t pi = 0; pi <= p; pi++)
        {
            if(pi == 0) Omg2 = Omg;
            if(std::pow(2, pi) >= batchs)
            {
                // reset H1, H2 to zero
                H1.setZero();
                H2.setZero();
            }
            band = std::fmin(band * 2, batchs);
            for(uint32_t b = 0, i = 1, j = 1; b < batchs; ++b, ++i, ++j)
            {
                start_idx = b * blocksize;
                stop_idx = (b + 1) * blocksize >= nrow ? nrow - 1 : (b + 1) * blocksize - 1;
                actual_block_size = stop_idx - start_idx + 1;
                if(trans)
                    G.middleRows(start_idx, actual_block_size).noalias() =
                        mat.middleCols(start_idx, actual_block_size).transpose() * Omg;
                else
                    G.middleRows(start_idx, actual_block_size).noalias() =
                        mat.middleRows(start_idx, actual_block_size) * Omg;
                if(i <= band / 2)
                {
                    if(trans)
                        H1.noalias() += mat.middleCols(start_idx, actual_block_size)
                                        * G.middleRows(start_idx, actual_block_size);
                    else
                        H1.noalias() += mat.middleRows(start_idx, actual_block_size).transpose()
                                        * G.middleRows(start_idx, actual_block_size);
                }
                else
                {
                    if(trans)
                        H2.noalias() += mat.middleCols(start_idx, actual_block_size)
                                        * G.middleRows(start_idx, actual_block_size);
                    else
                        H2.noalias() += mat.middleRows(start_idx, actual_block_size).transpose()
                                        * G.middleRows(start_idx, actual_block_size);
                }

                // use the first quarter band of succesive iteration (H1)
                // for extra power iteration updates with the last used band (H2)
                const bool adjacent = (pi > 0 && (b + 1) == std::pow(2, pi - 1) && std::pow(2, pi) < batchs);
                if((b + 1) < band && !adjacent) continue;
                if(!((i == band) || (i == band / 2) || adjacent)) continue;
                H = H1 + H2;
                Eigen::HouseholderQR<Matrix> qr(H);
                Omg.noalias() = qr.householderQ() * Matrix::Identity(cols(), size);
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

template<typename MatrixType, typename RsvdOp>
class RsvdOnePass
{
  private:
    using ConstGenericMatrix = const Eigen::Ref<const MatrixType>;
    using Index = Eigen::Index;
    using Scalar = typename MatrixType::Scalar;
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>; // dense matrix
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    Matrix b_leftSingularVectors{};
    Vector b_singularValues{};
    Matrix b_rightSingularVectors{};

  private:
    RsvdOp & b_op;

  public:
    RsvdOnePass(RsvdOp & op) : b_op(op) {}

    ~RsvdOnePass() {}

    inline Vector singularValues() const
    {
        return b_singularValues;
    }

    inline Matrix matrixU(bool trans = false) const
    {
        if(trans)
        {
            return b_rightSingularVectors;
        }
        else
        {
            return b_leftSingularVectors;
        }
    }

    inline Matrix matrixV(bool trans = false) const
    {
        if(trans)
        {
            return b_leftSingularVectors;
        }
        else
        {
            return b_rightSingularVectors;
        }
    }

    // G = D * Omega; H = D.transpose() * G;
    void computeUSV(uint32_t p, uint32_t batchs = 0)
    {
        const Eigen::Index nrow{b_op.rows()};
        const Eigen::Index ncol{b_op.cols()};
        const Eigen::Index size{b_op.ranks() + b_op.oversamples()};
        const Eigen::Index k{b_op.ranks()};
        Matrix H(ncol, size), G(nrow, size), R(size, size), Rt(size, size);

        if(batchs > 0)
            b_op.computeGandH(G, H, p, batchs);
        else
            b_op.computeGandH(G, H, p);

        {
            Eigen::HouseholderQR<Eigen::Ref<Matrix>> qr(G);
            R.noalias() =
                Matrix::Identity(size, nrow) * qr.matrixQR().template triangularView<Eigen::Upper>();
            G.noalias() = qr.householderQ() * Matrix::Identity(nrow, size);
        }
        {
            Eigen::HouseholderQR<Eigen::Ref<Matrix>> qr(G);
            Rt.noalias() =
                Matrix::Identity(size, nrow) * qr.matrixQR().template triangularView<Eigen::Upper>();
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
};

template<typename MatrixType>
class RsvdOne
{
  private:
    using ConstGenericMatrix = const Eigen::Ref<const MatrixType>;
    using Scalar = typename MatrixType::Scalar;
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>; // dense matrix
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;


    ConstGenericMatrix mat;
    const uint32_t k, os, rand;
    bool trans;

    RsvdOpOnePass<MatrixType> * op;
    RsvdOnePass<MatrixType, RsvdOpOnePass<MatrixType>> * rsvd;

  public:
    RsvdOne(ConstGenericMatrix & mat_, uint32_t k_, uint32_t os_ = 10, uint32_t rand_ = 1)
    : mat(mat_), k(k_), os(os_), rand(rand_)
    {
        if(mat.rows() >= mat.cols())
            trans = false;
        else
            trans = true;
        op = new RsvdOpOnePass<MatrixType>(mat, k, os, rand);
        rsvd = new RsvdOnePass<MatrixType, RsvdOpOnePass<MatrixType>>(*op);
    }

    ~RsvdOne()
    {
        delete op;
        delete rsvd;
    }

    inline void setRangeFinder(int flag)
    {
        op->finder = flag;
    }

    inline void compute(uint32_t p, uint32_t batchs = 0)
    {
        rsvd->computeUSV(p, batchs);
    }

    inline Matrix matrixU() const
    {
        return rsvd->matrixU(trans);
    }

    inline Matrix matrixV() const
    {
        return rsvd->matrixV(trans);
    }

    inline Vector singularValues() const
    {
        return rsvd->singularValues();
    }
};

template<typename MatrixType>
class RsvdDash
{
  private:
    using ConstGenericMatrix = const Eigen::Ref<const MatrixType>;
    using Index = Eigen::Index;
    using Scalar = typename MatrixType::Scalar;
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    ConstGenericMatrix mat;
    const uint32_t k, os, size, rand;
    const Index nrow, ncol;
    bool tall;
    Matrix Omg;
    Matrix b_leftSingularVectors;
    Vector b_singularValues;
    Matrix b_rightSingularVectors;

  public:
    RsvdDash(ConstGenericMatrix & mat_, uint32_t k_, uint32_t os_ = 10, uint32_t rand_ = 1)
    : mat(mat_), k(k_), os(os_), size(k_ + os_), rand(rand_), nrow(mat_.rows()), ncol(mat_.cols())
    {
        if(size > nrow || size > ncol)
            throw std::runtime_error(
                "k(pc) + s(oversampling) must be not greater than the dimension of data");
        tall = nrow > ncol ? true : false;
        auto randomEngine = std::default_random_engine{};
        if(rand == 1)
        {
            Omg = StandardNormalRandom<Matrix, std::default_random_engine>(tall ? nrow : ncol, size,
                                                                           randomEngine);
        }
        else
        {
            Omg = UniformRandom<Matrix, std::default_random_engine>(tall ? nrow : ncol, size, randomEngine);
        }
        b_leftSingularVectors = Matrix::Zero(nrow, size);
        b_rightSingularVectors = Matrix::Zero(ncol, size);
        b_singularValues = Vector::Zero(size);
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
        Omg = mat.transpose() * Omg;
        eigSVD(Omg, b_rightSingularVectors, b_singularValues, b_leftSingularVectors);
        double alpha = 0.0;
        for(uint32_t i = 0; i < p; i++)
        {
            Omg.noalias() = mat.transpose() * (mat * b_rightSingularVectors) - alpha * b_rightSingularVectors;
            eigSVD(Omg, b_rightSingularVectors, b_singularValues, b_leftSingularVectors);
            if(b_singularValues[size - 1] > alpha) alpha = (alpha + b_singularValues[size - 1]) / 2.0;
        }
        eigSVD(mat * b_rightSingularVectors, b_leftSingularVectors, b_singularValues, Omg);
        b_rightSingularVectors = b_rightSingularVectors * Omg;
    }

    void compute_wide(uint32_t p)
    {
        Omg = mat * Omg;
        eigSVD(Omg, b_leftSingularVectors, b_singularValues, b_rightSingularVectors);
        double alpha = 0.0;
        for(uint32_t i = 0; i < p; i++)
        {
            Omg.noalias() = mat * (mat.transpose() * b_leftSingularVectors) - alpha * b_leftSingularVectors;
            eigSVD(Omg, b_leftSingularVectors, b_singularValues, b_rightSingularVectors);
            if(b_singularValues[size - 1] > alpha) alpha = (alpha + b_singularValues[size - 1]) / 2.0;
        }
        eigSVD(mat.transpose() * b_leftSingularVectors, b_rightSingularVectors, b_singularValues, Omg);
        b_leftSingularVectors = b_leftSingularVectors * Omg;
    }

    /// A is tall and thin
    void eigSVD(const Matrix & A, Matrix & U, Vector & S, Matrix & V)
    {
        Matrix C = A.transpose() * A;
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
        const double tolerance = 1e-10 * S.maxCoeff() * A.cols();
        Vector inv_S = (S.array().abs() > tolerance).select(S.array().inverse(), 0.0);

        U.noalias() = A * V * inv_S.asDiagonal();
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
