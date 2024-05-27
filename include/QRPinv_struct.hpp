#ifndef ACE_QRPINV_STRUCT_DEFINED_H
#define ACE_QRPINV_STRUCT_DEFINED_H
#include <Eigen/Dense>

namespace ACE{

template <typename T> struct QRPinv_struct_T{
  //calculate A=Q*R*P as well as (R*P)^{-1}
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Q;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> RPinv;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> RPinv_inv;
  Eigen::VectorXd weights;

  void calculate(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & A);

  QRPinv_struct_T(){}
  QRPinv_struct_T(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & A){
    calculate(A);
  }
};

typedef QRPinv_struct_T<std::complex<double> > QRPinv_struct;
typedef QRPinv_struct_T<double > QRPinv_struct_real;

}//namespace
#endif
