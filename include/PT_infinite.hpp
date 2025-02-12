#ifndef ACE_PT_INFINITE_DEFINED_H
#define ACE_PT_INFINITE_DEFINED_H

#include <memory>
#include "ProcessTensorForward.hpp"
#include "Parameters.hpp"

namespace ACE{


struct infinite_normalize_specs{
  bool dont_normalize;
  bool use_iter;   //use any iterative method (Power iteration or Arnoldi)
  bool use_Arnoldi; 
  int iter; //maximum number of iterations
  double eps; //accuracy
};

Eigen::MatrixXcd PT_iTEBD_calc_R(const Eigen::VectorXcd & LambdaB, const Eigen::MatrixXcd &GLG2);

Eigen::MatrixXcd PT_iTEBD_calc_L(const Eigen::VectorXcd & LambdaB, const Eigen::MatrixXcd &GLG2);

Eigen::VectorXcd PT_iTEBD_R_times(const Eigen::VectorXcd & LambdaB, const Eigen::MatrixXcd &GLG2, const Eigen::VectorXcd &vr);

Eigen::VectorXcd PT_iTEBD_LT_times(const Eigen::VectorXcd & LambdaB, const Eigen::MatrixXcd &GLG2, const Eigen::VectorXcd &VL);

//PT_iTEBD_X_Result PT_iTEBD_X(const Eigen::VectorXcd & LambdaB, const Eigen::MatrixXcd &GLG2, const infinite_normalize_specs & specs);

//PT_iTEBD_X_Result PT_iTEBD_Y(const Eigen::VectorXcd & LambdaB, const Eigen::MatrixXcd &GLG2, const infinite_normalize_specs & specs);

void PT_iTEBD_step(Eigen::VectorXcd & LambdaA, Eigen::VectorXcd &LambdaB,
                   MPS_Matrix & GammaA, MPS_Matrix & GammaB,
                   Eigen::MatrixXcd expS, const TruncatedSVD &trunc, 
                   const infinite_normalize_specs & specs);

void PT_TEBD_step(Eigen::VectorXcd & LambdaA, const Eigen::VectorXcd &LambdaB,
                   MPS_Matrix & GammaA, MPS_Matrix & GammaB, 
                   const Eigen::MatrixXcd expS, const TruncatedSVD &trunc);


std::shared_ptr<ProcessTensorForward> PT_infinite(Parameters &param, DiagBB &diagBB);

}//namespace
#endif
