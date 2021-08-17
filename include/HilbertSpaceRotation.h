#ifndef ACE_HILBERT_SPACE_ROTATION_DEFINED_H
#define ACE_HILBERT_SPACE_ROTATION_DEFINED_H

#include "otimes.h"

/*
   Purpose: Store and apply rotations of the Hilbert space to propagators,
   initial values, and output operators
*/

class HilbertSpaceRotation{
public:

  bool use;  //do nothing, if false
  Eigen::MatrixXcd U; //Columns: Basis vectors

  bool used()const{return use;}

  Eigen::MatrixXcd apply(const Eigen::MatrixXcd & A)const{
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

  Eigen::MatrixXcd apply_Liouville(const Eigen::MatrixXcd & L)const{
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
  
  void setup_by_diagonalizing(const Eigen::MatrixXcd & A){
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(A);
    U=solver.eigenvectors();
    use=true;
  }
  
  HilbertSpaceRotation(){
    use=false;
  }

};



#endif
