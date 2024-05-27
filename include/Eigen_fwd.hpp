#pragma once
#ifndef ACE_EIGEN_FWD_DEFINED_H
#define ACE_EIGEN_FWD_DEFINED_H

#include <complex>

namespace Eigen {
    template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
    class Matrix;
    using MatrixXcd = Matrix<std::complex<double>, -1, -1, 0, -1, -1>;
    using VectorXcd = Matrix<std::complex<double>, -1, 1, 0, -1, 1>;
    using MatrixXd = Matrix<double, -1, -1, 0, -1, -1>;
    using VectorXd = Matrix<double, -1, 1, 0, -1, 1>;
    using Vector3d = Matrix<double, 3, 1, 0, 3, 1>;
}

#endif
