#ifndef MULTITIMEOP_DEFINED_H
#define MULTITIMEOP_DEFINED_H

#include "otimes.h"

class MultitimeOp{
public:
  double t;
  Eigen::MatrixXcd op;

  bool is_now(double t_, double dt){
    if( t_ >= t - 1e-8*dt  && t_ < t +dt - 1e-8*dt){
      return true;
    }else{
      return false;
    }
  }

  void apply(Eigen::MatrixXcd &M){
    if( op.rows() != op.cols() ){
      std::cerr<<"Applying Multitime Operator: op.rows() != op.cols() !"<<std::endl; 
      exit(1);
    } 
    if( M.cols() != op.rows() ){
      std::cerr<<"Applying Multitime Operator: M.cols() != op.rows() !"<<std::endl; 
      std::cerr<<"M.cols(): "<<M.cols()<<" op.rows(): "<<op.rows()<<std::endl;
      exit(1);
    } 
   
    M=M*op;
      
  }
    
  void apply(double t_, double dt, Eigen::MatrixXcd &M){
    if(is_now(t_,dt)){
      apply(M);
    }
  }

  void set(double t_, Eigen::MatrixXcd mat_fw, Eigen::MatrixXcd mat_bw){
    t=t_;
//    op=otimes(mat_fw, mat_bw.adjoint() );

    int dim=mat_fw.rows();
    Eigen::MatrixXcd id=Eigen::MatrixXcd::Identity(dim,dim);

    op=otimes(mat_fw, id) * otimes(id,mat_bw );
//    op=otimes(mat_fw, id) * otimes(id,mat_bw.adjoint() );

//std::cout<<"Test: Multitime: op: "<<std::endl<<op<<std::endl;
  } 
  void set(double t_, Eigen::MatrixXcd mat){
    t=t_;
    op=mat;
  } 
  MultitimeOp(){
  }
  MultitimeOp(double t_, Eigen::MatrixXcd mat_fw, Eigen::MatrixXcd mat_bw){ 
    set(t_, mat_fw, mat_bw);
  }
  MultitimeOp(double t_=0., Eigen::MatrixXcd mat=Eigen::MatrixXcd::Identity(2,2)){
    set(t_, mat);
  }
};


#endif
