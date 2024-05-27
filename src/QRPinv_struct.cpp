#include "QRPinv_struct.hpp"
#include <iostream>

namespace ACE{

template <typename T> 
void QRPinv_struct_T<T>::calculate(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & A){
//      Eigen::ColPivHouseholderQR<Eigen::MatrixXcd> qr(A);
    Eigen::FullPivHouseholderQR<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > qr(A);
    int rank=qr.rank();

    Q=qr.matrixQ().block(0,0,A.rows(),rank);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> R=
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(
          qr.matrixQR().template triangularView<Eigen::Upper>())
             .block(0,0,rank, A.cols() );

    Eigen::Matrix<T, Eigen::Dynamic, 1> v=R.diagonal().head(rank);
    weights=Eigen::VectorXd(rank); 
    for(int i=0; i<rank; i++){weights(i)=std::abs(v(i));}
//    weights=Eigen::VectorXd::Ones(rank);

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> P=qr.colsPermutation();
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Pinv=qr.colsPermutation().inverse();
    RPinv=R*Pinv;

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Rinv=
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(A.cols(),rank);
    Rinv.block(0,0,rank,rank)=qr.matrixQR().block(0,0,rank,rank).template triangularView<Eigen::Upper>().solve(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(rank,rank));

    RPinv_inv=P*Rinv;
  }

template class QRPinv_struct_T<std::complex<double> >;
template class QRPinv_struct_T<double>;

}//namespace
