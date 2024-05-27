#ifndef COMPRESS_TRAFO_AT_DEFINED_H
#define COMPRESS_TRAFO_AT_DEFINED_H

#include <Eigen/Core>
#include "MPS_Matrix.hpp"

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

namespace ACE{
class Compress_Trafo_At{
public:
  int n;
  Eigen::MatrixXcd R, L;
  bool use_ortho;
  
  void set_L_if_correct_n(const Eigen::MatrixXcd &L_, int n_);
  void set_R_if_correct_n(const Eigen::MatrixXcd &R_, int n_);
  
  Eigen::MatrixXcd get_Linv()const;
  Eigen::MatrixXcd get_Rinv()const;
  
  void apply_trafo(Eigen::MatrixXcd &M, bool low_high_low)const;
  void apply_trafo(MPS_Matrix &M, bool low_high_low)const;

  void apply_to_op(Eigen::VectorXcd &op, bool low_high_low)const;
  void apply_to_init(Eigen::VectorXcd &init, bool low_high_low)const;
  
  inline Compress_Trafo_At(int n_=-1) : n(n_){
    use_ortho=false;
  }
};

}//namespace
#endif
