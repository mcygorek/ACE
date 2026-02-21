#include "MultitimeOp.hpp"
#include "otimes.hpp"
#include "DummyException.hpp"
#include <iostream>

namespace ACE{

  bool MultitimeOp::is_now(double t_, double dt)const{
    if(apply_before){ 
      if( t_ >= (t-dt) - 1e-8*dt  && t_ < (t-dt) +dt - 1e-8*dt){
        return true;
      }else{
        return false;
      }
    }else{
      if( t_ >= t - 1e-8*dt  && t_ < t +dt - 1e-8*dt){
        return true;
      }else{
        return false;
      }
    }
  }

  void MultitimeOp::apply(Eigen::MatrixXcd &M){
    if( op_fw.rows() != op_fw.cols() ){
      std::cerr<<"Applying Multitime Operator: op_fw.rows() != op_fw.cols() !"<<std::endl; 
      throw DummyException();
    } 
    if( op_bw.rows() != op_bw.cols() ){
      std::cerr<<"Applying Multitime Operator: op_bw.rows() != op_bw.cols() !"<<std::endl; 
      throw DummyException();
    } 
    if( M.cols() != op_fw.rows()*op_bw.rows() ){
      std::cerr<<"Applying Multitime Operator: M.cols() != op.rows() !"<<std::endl; 
      std::cerr<<"M.cols(): "<<M.cols()<<" op.rows(): "<<op_fw.rows()*op_bw.rows()<<std::endl;
      throw DummyException();
    } 

    Eigen::MatrixXcd op=otimes(op_fw, op_bw.transpose() );

    if(apply_before){ 
      M=op*M;
    }else{
      M=M*op;
    }
  }
    
  void MultitimeOp::apply(double t_, double dt, Eigen::MatrixXcd &M){
    if(is_now(t_,dt)){
      apply(M);
    }
  }

  void MultitimeOp::set(double t_, Eigen::MatrixXcd mat_fw, Eigen::MatrixXcd mat_bw, bool before){
    t=t_;

    op_fw=mat_fw; 
    op_bw=mat_bw; 
    apply_before=before;

    if( op_fw.rows() != op_fw.cols() ){
      std::cerr<<"Setting Multitime Operator: op_fw.rows() != op_fw.cols() !"<<std::endl; 
      throw DummyException();
    } 
    if( op_bw.rows() != op_bw.cols() ){
      std::cerr<<"Setting Multitime Operator: op_bw.rows() != op_bw.cols() !"<<std::endl; 
      throw DummyException();
    } 
    if( op_fw.rows() != op_bw.rows() ){
      std::cerr<<"Settining Multitime Operator: op_fw.rows() != op_bw.rows() !"<<std::endl; 
      std::cerr<<"op_fw.rows(): "<<op_fw.rows()<<" op_bw.rows(): "<<op_bw.rows()*op_bw.rows()<<std::endl;
      throw DummyException();
    } 
  } 
   
}//namespace
