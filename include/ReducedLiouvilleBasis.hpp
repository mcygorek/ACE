#ifndef ACE_REDUCED_LIOUVILLE_BASIS_DEFINED_H
#define ACE_REDUCED_LIOUVILLE_BASIS_DEFINED_H

#include <Eigen/Dense>

/** Given a single multi-dimensional environment mode, there may be
more or less compact representation capturing the dynamics 
(e.g. coherent states for bosons) depending potentially on the initial
state of the environment. 
This interface can be used to link the original basis to the reduced basis.
We work in Liouville space with a matrix U containing relevant vectors as 
columns.
*/
namespace ACE{

class ReducedLiouvilleBasis{
public:

  bool use_reduce;
  Eigen::MatrixXcd U;   

  inline bool use()const{ return use_reduce; }
  inline virtual int get_reduced_dim()const{ return U.cols(); }
  inline virtual int get_original_dim()const{ return U.rows(); }

  inline int override_dim(int ML)const{
    if(use())return get_reduced_dim();
    return ML;
  }

  //reduce dimension of environment part of propagator: U^\dagger Prop_E U:
  Eigen::MatrixXcd transform_prop(const Eigen::MatrixXcd &mat, int Ns=0)const;

  //reduce dimension of initial 
  Eigen::VectorXcd transform_init(const Eigen::VectorXcd &vec)const;
  
  //reduce dimension of operator (e.g. env_op):
  Eigen::VectorXcd transform_op(const Eigen::VectorXcd &vec)const;
  
  //get reduced identity operator:
  Eigen::VectorXcd id_op(int M)const;

  inline ReducedLiouvilleBasis(const Eigen::MatrixXcd &U_): U(U_){
    use_reduce=true;
  }
  inline ReducedLiouvilleBasis(){
    use_reduce=false;
  }
  inline virtual ~ReducedLiouvilleBasis(){}
};
}//namespace
#endif
