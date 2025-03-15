/***********************************************************************
 * @file        Rand.hpp
 * @author      Zilong Li
 * @email       zilong.dk@gmail.com
 * @brief       Various Randomized Singular Value Decomposition in Eigen
 * Copyright (C) 2021 Zilong Li
 ***********************************************************************/

#ifndef UTIL_RAND_H
#define UTIL_RAND_H

#include <RcppEigen.h>
// #include <Eigen/Dense>
#include <random>

template<typename MatrixType>
using RealType = typename Eigen::NumTraits<typename MatrixType::Scalar>::Real;

template<typename MatrixType, typename ScalarType, typename RandomEngineType>
struct StandardNormalRandomHelper
{
    static inline MatrixType generate(Eigen::Index numRows, Eigen::Index numCols, RandomEngineType & engine);
};

template<typename MatrixType, typename RandomEngineType>
struct StandardNormalRandomHelper<MatrixType, RealType<MatrixType>, RandomEngineType>
{
    static inline MatrixType generate(const Eigen::Index numRows,
                                      const Eigen::Index numCols,
                                      RandomEngineType & engine)
    {
        std::normal_distribution<RealType<MatrixType>> distribution{0, 1};
        const auto normal{[&](typename MatrixType::Scalar) { return distribution(engine); }};
        return MatrixType::NullaryExpr(numRows, numCols, normal);
    }
};

template<typename MatrixType, typename RandomEngineType>
struct StandardNormalRandomHelper<MatrixType, std::complex<RealType<MatrixType>>, RandomEngineType>
{
    static inline MatrixType generate(const Eigen::Index numRows,
                                      const Eigen::Index numCols,
                                      RandomEngineType & engine)
    {
        constexpr RealType<MatrixType> stdDev{0.707106781186547};
        std::normal_distribution<RealType<MatrixType>> distribution{0, stdDev};
        const auto complexNormal{[&](typename MatrixType::Scalar) {
            return std::complex<RealType<MatrixType>>(distribution(engine), distribution(engine));
        }};
        return MatrixType::NullaryExpr(numRows, numCols, complexNormal);
    }
};

template<typename MatrixType, typename RandomEngineType>
inline MatrixType UniformRandom(const Eigen::Index numRows,
                                const Eigen::Index numCols,
                                RandomEngineType & engine)
{
    std::uniform_real_distribution<typename MatrixType::Scalar> uniform_real_distribution{-1, 1};
    const auto uniform{[&](typename MatrixType::Scalar) { return uniform_real_distribution(engine); }};
    return MatrixType::NullaryExpr(numRows, numCols, uniform);
}

template<typename MatrixType, typename RandomEngineType>
inline MatrixType StandardNormalRandom(const Eigen::Index numRows,
                                       const Eigen::Index numCols,
                                       RandomEngineType & engine)
{
    return StandardNormalRandomHelper<MatrixType, typename MatrixType::Scalar, RandomEngineType>::generate(
        numRows, numCols, engine);
}

#endif // UTIL_RAND_H
