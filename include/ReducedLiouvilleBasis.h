#ifndef ACE_REDUCED_LIOUVILLE_BASIS_DEFINED_H
#define ACE_REDUCED_LIOUVILLE_BASIS_DEFINED_H

/** Given a single multi-dimensional environment mode, there may be
more or less compact representation capturing the dynamics 
(e.g. coherent states for bosons) depending potentially on the initial
state of the environment. 
This interface can be used to link the original basis to the reduced basis.
We work in Liouville space with a matrix U containing relevant vectors as 
columns.
*/

class ReducedLiouvilleBasis{
public:

  bool use_reduce;
  Eigen::MatrixXcd U;   

  bool use()const{ return use_reduce; }
  virtual int get_reduced_dim()const{ return U.cols(); }
  virtual int get_original_dim()const{ return U.rows(); }

  int override_dim(int ML)const{
    if(use())return get_reduced_dim();
    return ML;
  }

  //reduce dimension of environment part of propagator: U^\dagger Prop_E U:
  Eigen::MatrixXcd transform_prop(const Eigen::MatrixXcd &mat, int Ns=0)const{
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
  Eigen::VectorXcd transform_init(const Eigen::VectorXcd &vec)const{
    if(!use())return vec;
    if(vec.size()!=U.rows()){
      std::cerr<<"ReducedLiouvilleBasis::transform_init: vec.size()!=U.rows()!"<<std::endl;
      exit(1);
    }
    return U.adjoint()*vec;
  }
  //reduce dimension of operator (e.g. env_op):
  Eigen::VectorXcd transform_op(const Eigen::VectorXcd &vec)const{
    if(!use())return vec;
    if(vec.size()!=U.rows()){
      std::cerr<<"ReducedLiouvilleBasis::transform_op: vec.size()!=U.rows()!"<<std::endl;
      exit(1);
    }
    return vec.transpose()*U;
  }
  //get reduced identity operator:
  Eigen::VectorXcd id_op(int M)const{
    Eigen::VectorXcd v=H_Matrix_to_L_Vector(Eigen::MatrixXcd::Identity(M,M));
    return transform_op(v);
  }


  ReducedLiouvilleBasis(const Eigen::MatrixXcd &U_): U(U_){
    use_reduce=true;
  }
  ReducedLiouvilleBasis(){
    use_reduce=false;
  }
  virtual ~ReducedLiouvilleBasis(){}
};



#endif
