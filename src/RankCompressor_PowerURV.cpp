#include "RankCompressor_PowerURV.hpp"
#include "PowerURV.hpp"

namespace ACE{

  void RankCompressor_PowerURV::compress(Eigen::MatrixXcd &A, Eigen::MatrixXcd &L, Eigen::MatrixXcd &R, bool left_to_right){

    PowerURV(A,L,R,threshold,debug);

  if(false){
//  if(true){
    if(left_to_right){
      if(A.rows()<A.cols()){  
        // A -> (A^dagger)^dagger -> (L' R')^dagger -> (R')^dagger (L')^dagger
        // move weights (contained in R') to right hand side
        for(int l=0; l<L.cols(); l++){
          double n=L.col(l).norm();
          if(n>1e-18){
            L.col(l)/=n;
            R.row(l)*=n;
          }
        }
      }
    }else{
      if(A.rows()>A.cols()){
        for(int l=0; l<L.cols(); l++){
          double n=R.col(l).norm();
          if(n>1e-18){
            R.row(l)/=n;
            L.col(l)*=n;
          }
        }
      }
    }
  }

/*
    std::cout<<"A.rows(): "<<A.rows()<<" A.cols(): "<<A.cols()<<" L.rows(): "<<L.rows()<<" rank: "<<L.cols()<<" R.cols(): "<<R.cols()<<std::endl;
    if(L.cols()<1){
      std::cout<<"A:"<<std::endl<<A<<std::endl<<std::endl;
      exit(1);
    }
*/
  }

}//namespace
