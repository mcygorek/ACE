#include "MultitimeOp.hpp"
#include "otimes.hpp"
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
    if( op.rows() != op.cols() ){
      std::cerr<<"Applying Multitime Operator: op.rows() != op.cols() !"<<std::endl; 
      exit(1);
    } 
    if( M.cols() != op.rows() ){
      std::cerr<<"Applying Multitime Operator: M.cols() != op.rows() !"<<std::endl; 
      std::cerr<<"M.cols(): "<<M.cols()<<" op.rows(): "<<op.rows()<<std::endl;
      exit(1);
    } 
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

    int dim=mat_fw.rows();
    Eigen::MatrixXcd id=Eigen::MatrixXcd::Identity(dim,dim);

    op=otimes(mat_fw, id) * otimes(id,mat_bw.transpose() );

    apply_before=before;
//std::cout<<"Test: Multitime: op: "<<std::endl<<op<<std::endl;
  } 
   
}//namespace
