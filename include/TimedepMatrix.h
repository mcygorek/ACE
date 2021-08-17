#ifndef ACE_TIMEDEP_MATRIX_DEFINED_H
#define ACE_TIMEDEP_MATRIX_DEFINED_H

#include <unsupported/Eigen/MatrixFunctions>
#include "Function.h"
#include "CheckMatrix.h"

class TimedepMatrix{
public:
  virtual Eigen::MatrixXcd f(double t) const =0;
  virtual int get_dim()const{
    return f(0).rows();
  }

  virtual ~TimedepMatrix(){}
};

typedef Smart_Ptr<TimedepMatrix> TimedepMatrixPtr;



class TimedepMatrix_Const: public TimedepMatrix{
public:
  Eigen::MatrixXcd A;
 
  virtual int get_dim()const{
    return A.rows();
  }
  virtual Eigen::MatrixXcd f(double t) const{
    return A;
  }

  TimedepMatrix_Const(const Eigen::MatrixXcd A_) :  A(A_){
  }
  TimedepMatrix_Const(){}
};

class TimedepMatrix_Const_Shape: public TimedepMatrix{
public:
  ComplexFunctionPtr fct;
  Eigen::MatrixXcd A;
 
  virtual int get_dim()const{
    return A.rows();
  }
  virtual Eigen::MatrixXcd f(double t) const{
    return fct->f(t)*A;
  }

  TimedepMatrix_Const_Shape(ComplexFunctionPtr &fct_, 
                            const Eigen::MatrixXcd A_) 
    :  fct(fct_), A(A_){
  }
  TimedepMatrix_Const_Shape(){}
};

class TimedepMatrix_Const_Shape_Hermitian: public TimedepMatrix{
public:
  ComplexFunctionPtr fct;
  Eigen::MatrixXcd A;
 
  virtual int get_dim() const{
    return A.rows();
  }
  virtual Eigen::MatrixXcd f(double t) const{
    Eigen::MatrixXcd B=fct->f(t)*A;
    return B+B.adjoint();
  }

  TimedepMatrix_Const_Shape_Hermitian(ComplexFunctionPtr &fct_, 
                            const Eigen::MatrixXcd A_) 
    :  fct(fct_), A(A_){
  }
  TimedepMatrix_Const_Shape_Hermitian(){}
};

class TimedepMatrix_InteractionPicture: public TimedepMatrix{
public:
  Eigen::MatrixXcd H0, A;
 
  virtual int get_dim() const{
    return A.rows();
  }
  virtual Eigen::MatrixXcd f(double t) const{
    Eigen::MatrixXcd B=std::complex<double>(0,t/Constants::hbar_in_meV_ps)*H0;
    Eigen::MatrixXcd C=B.exp();
    return C*A*(C.adjoint());
  }

  TimedepMatrix_InteractionPicture(const Eigen::MatrixXcd H0_,
                                   const Eigen::MatrixXcd A_) 
     :  H0(H0_), A(A_){

      int d1=get_dim_check_square(H0, "TimedepMatrix_InteractionPicture: H0");
      check_matrix_square(A, "TimedepMatrix_InteractionPicture: A");
      check_matrix_rows_eq(A, d1, "TimedepMatrix_InteractionPicture: A");
     
  }
  TimedepMatrix_InteractionPicture(){}
};




#endif
