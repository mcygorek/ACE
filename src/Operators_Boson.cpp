#include "PCH.hpp"
#include <Eigen/Core>
#include "Operators_Boson.hpp"

namespace ACE{
namespace Operators_Boson{

  Eigen::MatrixXcd a(int n_max){
    Eigen::MatrixXcd result=Eigen::MatrixXcd::Zero(n_max, n_max);
    for(int i=0; i<n_max-1; i++){
      result(i, i+1)=sqrt(i+1);
    }
    return result;
  }

  Eigen::MatrixXcd adagger(int n_max){
    Eigen::MatrixXcd result=Eigen::MatrixXcd::Zero(n_max, n_max);
    for(int i=0; i<n_max-1; i++){
      result(i+1, i)=sqrt(i+1);
    }
    return result;
  }

  Eigen::MatrixXcd n(int n_max){
    Eigen::MatrixXcd result=Eigen::MatrixXcd::Zero(n_max, n_max);
    for(int i=0; i<n_max; i++){
      result(i, i)=i;
    }
    return result;
  } 
  Eigen::MatrixXcd id(int N){ 
    return Eigen::MatrixXcd::Identity(N,N); 
  }

  Eigen::MatrixXcd vacuum(int n_max){
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(n_max, n_max);
    mat(0,0)=1;
    return mat;
  }
};

//
namespace Operators_Boson_Offset{

  Eigen::MatrixXcd a(int base, int M){
    Eigen::MatrixXcd result=Eigen::MatrixXcd::Zero(M, M);
    for(int i=0; i<M-1; i++){
      result(i, i+1)=sqrt(base+i+1);
    }
    return result;
  }

  Eigen::MatrixXcd adagger(int base, int M){
    Eigen::MatrixXcd result=Eigen::MatrixXcd::Zero(M,M);
    for(int i=0; i<M-1; i++){
      result(i+1, i)=sqrt(base+i+1);
    }
    return result;
  }

  Eigen::MatrixXcd n(int base, int M){
    Eigen::MatrixXcd result=Eigen::MatrixXcd::Zero(M, M);
    for(int i=0; i<M; i++){
      result(i, i)=base+i;
    }
    return result;
  } 

  Eigen::MatrixXcd lowest(int M){
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(M, M);
    mat(0,0)=1;
    return mat;
  }
};

}//namespace
