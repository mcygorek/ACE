#include "ApproxSVD.hpp"
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/src/misc/RealSvd2x2.h>
#include <cstdlib>
#include <vector>
#include <algorithm>

namespace ACE{

template <typename ScalarType> double conditionalReal(ScalarType z){
  return z;
}
template<> double conditionalReal<std::complex<double> >(std::complex<double> z){
  return z.real();
}


template<typename ScalarType> 
  void ApproxSVD<ScalarType>::rotate_single(int p, int q){

    Eigen::JacobiRotation<double> j_left,j_right;
//    Eigen::internal::real_2x2_jacobi_svd<Eigen::MatrixType, double, int>(M, p, q, &j_left, &j_right);
    Eigen::internal::real_2x2_jacobi_svd(M, p, q, &j_left, &j_right);

    M.applyOnTheLeft(p,q,j_left);
    L.applyOnTheRight(p,q,j_left.transpose());

    M.applyOnTheRight(p,q,j_right);
    R.applyOnTheLeft(p,q,j_right.transpose());
  }

template<typename ScalarType> 
  void ApproxSVD<ScalarType>::initialize(const MatrixType &A, bool forceQR){
    if(!forceQR && A.rows()==A.cols()){
      M=A; 
      L=MatrixType::Identity(M.rows(), M.rows());
      R=MatrixType::Identity(M.cols(), M.cols());
    }else if(A.rows()>=A.cols()){
      int mbig=A.rows(); int msmall=A.cols();
      Eigen::ColPivHouseholderQR<MatrixType> qr;
      qr.compute(A);
      M=((MatrixType)qr.matrixR().template triangularView<Eigen::Upper>()).block(0,0,msmall,msmall); 
      L=((MatrixType)qr.householderQ()).block(0, 0, mbig, msmall);
      QRremainder=((MatrixType)qr.householderQ()).block(0, msmall, mbig, mbig-msmall);
      R=((MatrixType)qr.colsPermutation().transpose()).block(0,0,msmall, msmall);
    }else{
      int mbig=A.cols(); int msmall=A.rows();
      Eigen::ColPivHouseholderQR<MatrixType> qr;
      qr.compute(A.transpose());
      M=((MatrixType)qr.matrixR().template triangularView<Eigen::Upper>()).block(0,0,msmall,msmall).transpose(); 
      R=((MatrixType)qr.householderQ()).block(0, 0, mbig, msmall).transpose();
      QRremainder=((MatrixType)qr.householderQ()).block(0, msmall, mbig, mbig-msmall).transpose();
      L=((MatrixType)qr.colsPermutation()).block(0,0,msmall, msmall);
//      L=((MatrixType)qr.colsPermutation().transpose()).block(0,0,msmall, msmall).transpose();
    }
  }

  //In case L or R have been truncated by QR: reconstruct full matrices
template<typename ScalarType>
  typename ApproxSVD<ScalarType>::MatrixType ApproxSVD<ScalarType>::reconstructL(){  
    if(L.rows()==L.cols())return L;
    MatrixType L2(L.rows(),L.rows());
    L2.block(0,0,L.rows(),L.cols())=L;
    L2.block(0,L.cols(), QRremainder.rows(), QRremainder.cols())=QRremainder;
    return L2;
  }
template<typename ScalarType>
  typename ApproxSVD<ScalarType>::MatrixType ApproxSVD<ScalarType>::reconstructR(){    
    if(R.rows()==R.cols())return R;
    MatrixType R2(R.cols(),R.cols());
    R2.block(0,0,R.rows(),R.cols())=R;
    R2.block(R.rows(), 0, QRremainder.rows(), QRremainder.cols())=QRremainder;
    return R2;
  }

template<typename ScalarType> 
  bool ApproxSVD<ScalarType>::sort_second_descending(const std::pair<int, double> & first, 
                                     const std::pair<int, double> & second){
    return first.second > second.second;
  } 

template<typename ScalarType> 
  void ApproxSVD<ScalarType>::sort(){
    std::vector<std::pair<int, double> > vec(M.rows());
    for(int c=0; c<M.rows(); c++){
      vec[c]=std::make_pair(c,conditionalReal(M(c,c)));
    }
    std::sort(vec.begin(), vec.end(), sort_second_descending);

    {
    MatrixType L2(L.rows(), L.cols()); 
    for(int c=0; c<vec.size(); c++){
      L2.col(c)=L.col(vec[c].first);
    }
    L.swap(L2);
    }

    {
    MatrixType R2(R.rows(), R.cols()); 
    for(int c=0; c<vec.size(); c++){
      R2.row(c)=R.row(vec[c].first);
    }
    R.swap(R2);
    }

    {
    MatrixType M2(M.rows(), M.cols()); 
    for(int c=0; c<vec.size(); c++){
      M2.col(c)=M.col(vec[c].first);
    }
    for(int c=0; c<vec.size(); c++){
      M.row(c)=M2.row(vec[c].first);
    }
    }
  }

template<typename ScalarType> 
  void ApproxSVD<ScalarType>::calculate(double eps, int max_loop){
    if(M.rows()<2||M.cols()<2)return; //nothing to be done
    if(max_loop<1)max_loop=M.rows()*M.cols()*20;

    //Make a list of largest modulus off-diagonal element in each column
    Eigen::VectorXi max_in_c=Eigen::VectorXi::Zero(M.cols());   
    max_in_c(0)=1;
    for(int c=0; c<M.cols(); c++){
      for(int r=1; r<M.rows(); r++){
        if(r==c)continue;
        if(abs(M(r,c))>abs(M(max_in_c(c),c)))max_in_c(c)=r;
      }
    }


    //main loop:
    bool done=false; 
    int counter=0;
    while(!done){
#ifdef DEBUG_APPROXSVD
std::cout<<"loop: "<<counter<<std::endl;
#endif
      //find column with maximal maximal element
      int c_max=0;
      for(int c=1; c<M.cols(); c++){
        if(abs(M(max_in_c(c),c))>abs(M(max_in_c(c_max), c_max)))c_max=c;
      }
#ifdef DEBUG_APPROXSVD
std::cout<<"maximal element at ("<<max_in_c(c_max)<<", "<<c_max<<"): "<<M(max_in_c(c_max), c_max)<<std::endl;
#endif
      
      if(abs(M(max_in_c(c_max), c_max))<eps){
        done=true;
        break;
      }
    
      //eliminate maximal off-diagonal element
      int rot_r=max_in_c(c_max);
      int rot_c=c_max;
      rotate_single(rot_r, rot_c);

      //update list of largest modulus off-diagonal element in each column:
      //modified columns rot_r and rot_c: search brute force for new max:
      max_in_c(rot_r)= rot_r==0 ? 1 : 0;
      for(int r=0; r<M.rows(); r++){
        if(r==rot_r)continue;
        if(abs(M(r,rot_r))>abs(M(max_in_c(rot_r),rot_r)))max_in_c(rot_r)=r;
      }
      max_in_c(rot_c)= rot_c==0 ? 1 : 0;
      for(int r=0; r<M.rows(); r++){
        if(r==rot_c)continue;
        if(abs(M(r,rot_c))>abs(M(max_in_c(rot_c),rot_c)))max_in_c(rot_c)=r;
      }
      //for all others:
      for(int c=0; c<M.cols(); c++){
        if(c==rot_r||c==rot_c)continue;  //cases treated earlier
        if(max_in_c(c)!=rot_r && max_in_c(c)!=rot_c)continue; //still max.
        //redo from scratch:
        max_in_c(c)= c==0 ? 1 : 0;
        for(int r=0; r<M.rows(); r++){
          if(r==c)continue;
          if(abs(M(r,c))>abs(M(max_in_c(c),c)))max_in_c(c)=r;
        }
      }

      //emergency break to avoid infinite loops
      counter++;
      if(counter>max_loop){
        std::cout<<"WARNING: ApproxSVD has reached max_loop="<<max_loop<<"!"<<std::endl;
        break;
      }
    }
    // make sign of diagonals of M positive
    for(int c=0; c<M.cols(); c++){
      if(conditionalReal(M(c,c))<0){
        M.row(c)*=(-1.);
        L.col(c)*=(-1.);
      }
    }
    sort();
  }

template<typename ScalarType> 
  Eigen::VectorXd ApproxSVD<ScalarType>::singularValues()const{
    return M.diagonal().real();
  }

template<typename ScalarType> 
  void ApproxSVD<ScalarType>::swap(ApproxSVD &other){
    L.swap(other.L);
    M.swap(other.M);
    R.swap(other.R);
    QRremainder.swap(other.QRremainder);
  }

template<typename ScalarType> 
  void ApproxSVD<ScalarType>::compute(const MatrixType &A, double eps, bool forceQR){
    initialize(A, forceQR);
    calculate(eps);
  }


template<> void ApproxSVD<std::complex<double> >::rotate_single(int p, int q){
  typedef std::complex<double> Scalar;
  typedef double RealScalar;
  Scalar z;
  Eigen::JacobiRotation<Scalar> rot;
  RealScalar n = std::sqrt(Eigen::numext::abs2(M.coeff(p,p)) + Eigen::numext::abs2(M.coeff(q,p)));

  const RealScalar considerAsZero = (std::numeric_limits<RealScalar>::min)();

  if(n==0){
    // make sure first column is zero
    M.coeffRef(p,p) = M.coeffRef(q,p) = Scalar(0);

    if(std::abs(Eigen::numext::imag(M.coeff(p,q)))>considerAsZero){
      // M.coeff(p,q) can be zero if M.coeff(q,p) is not zero but small enough to underflow when computing n
      z = std::abs(M.coeff(p,q)) / M.coeff(p,q);
      M.row(p) *= z;
      L.col(p) *= std::conj(z);
    }
    if(abs(Eigen::numext::imag(M.coeff(q,q)))>considerAsZero){
      z = std::abs(M.coeff(q,q)) / M.coeff(q,q);
      M.row(q) *= z;
      L.col(q) *= std::conj(z);
    }
  // otherwise the second row is already zero, so we have nothing to do.
  }else{
    rot.c() = std::conj(M.coeff(p,p)) / n;
    rot.s() = M.coeff(q,p) / n;
    M.applyOnTheLeft(p,q,rot);
    L.applyOnTheRight(p,q,rot.adjoint());
    if(std::abs(Eigen::numext::imag(M.coeff(p,q)))>considerAsZero){
      z = std::abs(M.coeff(p,q)) / M.coeff(p,q);
      M.col(q) *= z;
      R.row(q) *= std::conj(z);
//      if(svd.computeV()) svd.m_matrixV.col(q) *= z;
    }
    if(std::abs(Eigen::numext::imag(M.coeff(q,q)))>considerAsZero){
      z = std::abs(M.coeff(q,q)) / M.coeff(q,q);
      M.row(q) *= z;
      L.col(q) *= std::conj(z);
    }
  }

  // update largest diagonal entry
//  maxDiagEntry = Eigen::numext::maxi<RealScalar>(maxDiagEntry,Eigen::numext::maxi<RealScalar>(std::abs(M.coeff(p,p)), std::abs(M.coeff(q,q))));
  // and check whether the 2x2 block is already diagonal
//  RealScalar threshold = Eigen::numext::maxi<RealScalar>(considerAsZero, precision * maxDiagEntry);

//  if(std::abs(M.coeff(p,q))>threshold || std::abs(M.coeff(q,p))>threshold){
  if(true){
    Eigen::JacobiRotation<double> j_left,j_right;
    Eigen::internal::real_2x2_jacobi_svd
        <Eigen::MatrixXcd,  double, int>  (M, p, q, &j_left, &j_right);

    M.applyOnTheLeft(p,q,j_left);
    L.applyOnTheRight(p,q,j_left.transpose());

    M.applyOnTheRight(p,q,j_right);
    R.applyOnTheLeft(p,q,j_right.transpose());
  }
}

template class ApproxSVD<double>;
template class ApproxSVD<std::complex<double>>;

}//namespace
