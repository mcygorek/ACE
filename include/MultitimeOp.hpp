#ifndef MULTITIMEOP_DEFINED_H
#define MULTITIMEOP_DEFINED_H

#include <Eigen/Dense>

namespace ACE{

class MultitimeOp{
public:
  double t;
  Eigen::MatrixXcd op;
  bool apply_before; //apply als last action of previous time step

  bool is_now(double t_, double dt)const;

  void apply(Eigen::MatrixXcd &M);
    
  void apply(double t_, double dt, Eigen::MatrixXcd &M);

  void set(double t_, Eigen::MatrixXcd mat_fw, Eigen::MatrixXcd mat_bw, bool before=false);
   
  inline void set(double t_, Eigen::MatrixXcd mat, bool before=false){
    t=t_;
    op=mat;
    apply_before=before;
  } 
  
  inline MultitimeOp(){
  }
  inline MultitimeOp(double t_, Eigen::MatrixXcd mat_fw, Eigen::MatrixXcd mat_bw, bool before=false){ 
    set(t_, mat_fw, mat_bw, before);
  }
  inline MultitimeOp(double t_=0., Eigen::MatrixXcd mat=Eigen::MatrixXcd::Identity(2,2), bool before=false){
    set(t_, mat, before);
  }
};

}//namespace
#endif
