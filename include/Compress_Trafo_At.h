#ifndef COMPRESS_TRAFO_AT_DEFINED_H
#define COMPRESS_TRAFO_AT_DEFINED_H

#include <Eigen/Core>


/** Sometimes, when compressing MPOs, we want to track what kind of 
    transformation has been used to compress a given chain element.
    This structure stors the respective transformation matrices.  

    E.g.: n=1:

    a[0] a[1] a[2] a[3] ... 
              ^     
              Track trafo between at this point
forward:
->  L[0] R[0]a[1] a[2] a[3] ...
->  L[0] L[1] R[1]a[2] a[3] ...   (=> save R[1])
->  L[0] L[1] L[2] R[2]a[3] ...   
... b[0] b[1] b[2] b[3]  ...
backward:
->  b[0] b[1] b[2]L[3] R[3] ...   
->  b[0] b[1]L[2] R[2] R[3] ...   (=> save L[2])
... c[0] c[1] c[2] c[3]  ...


Then, an operator O in the basis of the a's needs to be transformed:

    a[0] a[1] O a[2] a[3] ... 
->  L[0] L[1] R[1]*O*R^{-1}[1] R[1]a[2] a[3] ...
...
->  b[0] b[1] O' b[2]L[3] R[3] 
..
->  b[0] b[1]L[2] L^{-1}[2]*O'*L[2] R[2] R[3] ...       
... c[0] c[1] O'' c[2] c[3] ..

=> O'' = L^{-1}[n+1]*R[n]*O*R^{-1}[n]*L[n+1] 

(or, if sweels are from high to low to high
=> O'' = R[n]*L^{-1}[n+1]*O*L[n+1]*R^{-1}[n]
*/

class Compress_Trafo_At{
public:
  int n;
  Eigen::MatrixXcd R, L;
  bool use_ortho;
  
  void set_L_if_correct_n(const Eigen::MatrixXcd &L_, int n_){
    if(n_==n+1){
      L=L_;
      
      if(use_ortho)for(int c=0; c<L.cols(); c++){
        L.col(c).normalize();
      }
    }
  }
  void set_R_if_correct_n(const Eigen::MatrixXcd &R_, int n_){
    if(n_==n){
      R=R_;

      if(use_ortho)for(int r=0; r<R.rows(); r++){
        R.row(r).normalize();
      }
    }
  }
  Eigen::MatrixXcd get_Linv()const{
    Eigen::MatrixXcd Linv=L.adjoint();
    for(int r=0; r<L.cols(); r++){
      double norm=L.col(r).norm();
      Linv.row(r)/=norm*norm;
    }
    return Linv;
  }
  Eigen::MatrixXcd get_Rinv()const{
    Eigen::MatrixXcd Rinv=R.adjoint();
    for(int r=0; r<R.rows(); r++){
      double norm=R.row(r).norm();
      Rinv.col(r)/=norm*norm;
    }
    return Rinv;
  }
  void apply_trafo(Eigen::MatrixXcd &M, bool low_high_low)const{
    if(low_high_low){ 
      M=get_Linv()*R*M*get_Rinv()*L;
    }else{
      M=R*get_Linv()*M*L*get_Rinv();
    }
  }
  void apply_trafo(MPS_Matrix &M, bool low_high_low)const{
    if(low_high_low){
      M.inner_multiply_left(get_Linv()*R);
      M.inner_multiply_right(get_Rinv()*L); //M = Linv*R*M*Rinv*L;
    }else{
      M.inner_multiply_left(R*get_Linv());
      M.inner_multiply_right(L*get_Rinv()); //M = R*Linv*M*L*Rinv;
    }
  }
  void apply_to_op(Eigen::VectorXcd &op, bool low_high_low)const{
    if(low_high_low){
      op=get_Linv()*R*op;
    }else{
      std::cerr<<"Compress_Trafo_At::apply_to_op: not implemented!"<<std::endl;
      exit(1);
    }
  }
  void apply_to_init(Eigen::VectorXcd &init, bool low_high_low)const{
    if(low_high_low){
      init=(init.transpose()*get_Rinv()*L).transpose();
    }else{
      std::cerr<<"Compress_Trafo_At::apply_to_init: not implemented!"<<std::endl;
      exit(1);
    }
  }
  Compress_Trafo_At(int n_=-1) : n(n_){
    use_ortho=false;
  }
};

#endif
