#include "Compress_Trafo_At.hpp"
#include <iostream>
#include <Eigen/Core>

namespace ACE{
  
  void Compress_Trafo_At::set_L_if_correct_n(const Eigen::MatrixXcd &L_, int n_){
    if(n_==n+1){
      L=L_;
      
      if(use_ortho)for(int c=0; c<L.cols(); c++){
        L.col(c).normalize();
      }
    }
  }

  void Compress_Trafo_At::set_R_if_correct_n(const Eigen::MatrixXcd &R_, int n_){
    if(n_==n){
      R=R_;

      if(use_ortho)for(int r=0; r<R.rows(); r++){
        R.row(r).normalize();
      }
    }
  }

  Eigen::MatrixXcd Compress_Trafo_At::get_Linv()const{
    Eigen::MatrixXcd Linv=L.adjoint();
    for(int r=0; r<L.cols(); r++){
      double norm=L.col(r).norm();
      Linv.row(r)/=norm*norm;
    }
    return Linv;
  }

  Eigen::MatrixXcd Compress_Trafo_At::get_Rinv()const{
    Eigen::MatrixXcd Rinv=R.adjoint();
    for(int r=0; r<R.rows(); r++){
      double norm=R.row(r).norm();
      Rinv.col(r)/=norm*norm;
    }
    return Rinv;
  }

  void Compress_Trafo_At::apply_trafo(Eigen::MatrixXcd &M, bool low_high_low)const{
    if(low_high_low){ 
      M=get_Linv()*R*M*get_Rinv()*L;
    }else{
      M=R*get_Linv()*M*L*get_Rinv();
    }
  }

  void Compress_Trafo_At::apply_trafo(MPS_Matrix &M, bool low_high_low)const{
    if(low_high_low){
      M.inner_multiply_left(get_Linv()*R);
      M.inner_multiply_right(get_Rinv()*L); //M = Linv*R*M*Rinv*L;
    }else{
      M.inner_multiply_left(R*get_Linv());
      M.inner_multiply_right(L*get_Rinv()); //M = R*Linv*M*L*Rinv;
    }
  }

  void Compress_Trafo_At::apply_to_op(Eigen::VectorXcd &op, bool low_high_low)const{
    if(low_high_low){
      op=get_Linv()*R*op;
    }else{
      std::cerr<<"Compress_Trafo_At::apply_to_op: not implemented!"<<std::endl;
      exit(1);
    }
  }

  void Compress_Trafo_At::apply_to_init(Eigen::VectorXcd &init, bool low_high_low)const{
    if(low_high_low){
      Eigen::MatrixXcd Rinv=get_Rinv();
      if(init.rows()!=Rinv.rows()){
        std::cerr<<"Compress_Trafo_At::apply_to_init: init.rows()!=Rinv.rows() ("<<init.rows()<<"/"<<Rinv.rows()<<")!"<<std::endl;
        exit(1);
      }
      init=(init.transpose()*Rinv*L).transpose();
    }else{
      std::cerr<<"Compress_Trafo_At::apply_to_init: not implemented!"<<std::endl;
      exit(1);
    }
  }


}//namespace
