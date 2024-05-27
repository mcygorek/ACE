#include "PCH.hpp"
#include "HilbertSpaceRotation.hpp"
#include "otimes.hpp"
#include <iostream>
#include <cstdlib>


namespace ACE{

  Eigen::MatrixXcd HilbertSpaceRotation::apply(const Eigen::MatrixXcd & A)const{
    if(!used())return A;
    if(U.rows()!=A.cols()){
      std::cerr<<"HilbertSpaceRotation: U.rows()!=A.cols()!"<<std::endl;
      exit(1);
    }
    if(U.cols()!=A.rows()){
      std::cerr<<"HilbertSpaceRotation: U.cols()!=A.rows()!"<<std::endl;
      exit(1);
    }
    return U.adjoint()*A*U;
  }

  Eigen::MatrixXcd HilbertSpaceRotation::apply_Liouville(const Eigen::MatrixXcd & L)const{
    if(!used())return L;
    int N=U.rows();
    if(N*N!=L.cols()){
      std::cerr<<"HilbertSpaceRotation: N*N!=L.cols()!"<<std::endl;
      exit(1);
    }
    if(N*N!=L.rows()){
      std::cerr<<"HilbertSpaceRotation: N*N!=L.rows()!"<<std::endl;
      exit(1);
    }
    
    return otimes(U.adjoint(), U.transpose())* L *
           otimes(U, U.adjoint().transpose());
  }
  
  void HilbertSpaceRotation::setup_by_diagonalizing(const Eigen::MatrixXcd & A){
    bool is_diagonal=true;
    for(int r=0; r<A.rows(); r++){    
      for(int c=0; c<A.cols(); c++){    
        if(c!=r && abs(A(r,c))>1e-16){is_diagonal=false;}
      }
    }
    if(is_diagonal){
      use=false;
      U=Eigen::MatrixXcd::Identity(A.rows(), A.cols());
    }else{
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(A);
      U=solver.eigenvectors();
      use=true;
    }
  }

}//namespace
