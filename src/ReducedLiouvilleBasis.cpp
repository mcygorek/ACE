#include "ReducedLiouvilleBasis.hpp"
#include "LiouvilleTools.hpp"
#include "otimes.hpp"
#include <iostream> 

namespace ACE{

  //reduce dimension of environment part of propagator: U^\dagger Prop_E U:
  Eigen::MatrixXcd ReducedLiouvilleBasis::transform_prop(const Eigen::MatrixXcd &mat, int Ns)const{
    if(!use())return mat;
    if(Ns<2){
      if( mat.rows()!=U.rows() ){
        std::cerr<<"ReducedLiouvilleBasis::transform_prop: mat.rows()!=U.rows() ("<<mat.rows()<<" vs. "<<U.rows()<<")!"<<std::endl;
        exit(1);
      }
      return U.adjoint()*mat*U;
    }else{
      int NL=Ns*Ns;
      if( mat.rows()!=U.rows()*NL ){
        std::cerr<<"ReducedLiouvilleBasis::transform_prop: mat.rows()!=U.rows()*NL ("<<mat.rows()<<" vs. "<<U.rows()<<"*"<<NL<<")!"<<std::endl;
        exit(1);
      }
      Eigen::MatrixXcd V=otimes(Eigen::MatrixXcd::Identity(NL,NL), U);
      return V.adjoint()*mat*V;
    } 
  }

  //reduce dimension of initial 
  Eigen::VectorXcd ReducedLiouvilleBasis::transform_init(const Eigen::VectorXcd &vec)const{
    if(!use())return vec;
    if(vec.size()!=U.rows()){
      std::cerr<<"ReducedLiouvilleBasis::transform_init: vec.size()!=U.rows()!"<<std::endl;
      exit(1);
    }
    return U.adjoint()*vec;
  }

  //reduce dimension of operator (e.g. env_op):
  Eigen::VectorXcd ReducedLiouvilleBasis::transform_op(const Eigen::VectorXcd &vec)const{
    if(!use())return vec;
    if(vec.size()!=U.rows()){
      std::cerr<<"ReducedLiouvilleBasis::transform_op: vec.size()!=U.rows()!"<<std::endl;
      exit(1);
    }
    return vec.transpose()*U;
  }

  //get reduced identity operator:
  Eigen::VectorXcd ReducedLiouvilleBasis::id_op(int M)const{
    Eigen::VectorXcd v=H_Matrix_to_L_Vector(Eigen::MatrixXcd::Identity(M,M));
    return transform_op(v);
  }

}//namespace
