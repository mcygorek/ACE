#include "HermitianLiouvilleBasis.hpp"
#include <Eigen/Dense>

namespace ACE{

  Eigen::VectorXcd HermitianLiouvilleBasis::operator() (int i) const{

    if(i==0){  //proportional to trace
      Eigen::VectorXcd res=Eigen::VectorXcd::Zero(N*N);
      for(size_t j=0; j<N; j++){
        res(j*N+j)=1./sqrt((double) N);
      }
      return res;
    }else if(i<N){ //traceless diagonal (using DCT-II);
      Eigen::VectorXcd res=Eigen::VectorXcd::Zero(N*N);
      for(int n=0; n<N; n++){
        res(n*N+n)=sqrt(2./N)*cos(M_PI/N*(n+0.5)*i);
      }
      return res;
    }else if(i<N+((N-1)*N)/2){ 
      int j;
      int l=i-N;
      for(j=0; j<N-1; j++){
        if(l<(N-1-j)){
          l+=j+1; 
          break;
        }else{
          l-=(N-1-j);
        }
      }
      Eigen::VectorXcd res=Eigen::VectorXcd::Zero(N*N);
      res(j*N+l)=1./sqrt(2.);
      res(l*N+j)=1./sqrt(2.);
      return res;
    }else if(i<N*N){ 
      int j;
      int l=i-(N+((N-1)*N)/2);
      for(j=0; j<N-1; j++){
        if(l<(N-1-j)){
          l+=j+1; 
          break;
        }else{
          l-=(N-1-j);
        }
      }
      Eigen::VectorXcd res=Eigen::VectorXcd::Zero(N*N);
      res(j*N+l)=std::complex<double>(0., 1./sqrt(2.));
      res(l*N+j)=std::complex<double>(0.,-1./sqrt(2.));
      return res;
    }else{
      throw "HermitianLiouvilleBasis: out of bounds!";
    }
  }
  
  Eigen::MatrixXcd HermitianLiouvilleBasis::get_Matrix()const{
    Eigen::MatrixXcd M(N*N, N*N);
    for(int i=0; i<N*N; i++){
      M.col(i)=operator()(i);
    }
    return M;
  }

}//namespace
