#ifndef MPO_OP_DEFINED_H_
#define MPO_OP_DEFINED_H_

#include "MPO.hpp"

namespace ACE{
namespace MPO_op{  //features with clear interpretation of MPO as operator

  std::complex<double> evaluate(const MPO & mpo, const std::vector<int> &list_d1, const std::vector<int> &list_d3);

  std::complex<double> evaluate(const MPO & mpo, const std::vector<Eigen::VectorXcd> &closure_d1, const std::vector<int> &list_d3);

  std::vector<Eigen::VectorXcd> get_cap_d2(const MPO & mpo, const std::vector<Eigen::VectorXcd> &closure_d1);
  std::vector<Eigen::VectorXcd> get_cap_d4(const MPO & mpo, const std::vector<Eigen::VectorXcd> &closure_d1);

  void join_d3_d1(MPO & mpo, const MPO & other);
  void join_d3_d1_compress_d2_d4(MPO & mpo, const MPO & other, const TruncatedSVD &trunc, int verbosity=0);


  void add_coeff(MPO & mpo, const std::vector<int> &list_d1, const std::vector<int> &list_d3, const std::complex<double> &value);
 
  MPO H_to_L_forward(const MPO &mpo);
  MPO H_to_L_backward(const MPO &mpo);
  FourLeg join_H_to_L(const FourLeg & f1, const FourLeg & f2, bool conjugate_first=true);
  MPO join_H_to_L(const MPO & mpo1, const MPO & mpo2, bool conjugate_first=true);

  MPO SingleExcitationManifold_to_MPO_H(const Eigen::MatrixXcd & H, const TruncatedSVD & trunc, std::vector<int> dims=std::vector<int>());

  std::complex<double> SingleExcitationManifold_evaluate_MPO(const MPO & mpo, int i, int j);

  double compare_SingeExcitationManifold_MPO(const MPO &mpo, const Eigen::MatrixXcd &M);


  MPO Identity(const std::vector<int> & dims2);
  MPO Zero(const std::vector<int> & dims2);

  MPO OneBody(const std::vector<int> & dims, int r1, const Eigen::MatrixXcd &op1);
  MPO OneBody_forward(const std::vector<int> & dims2, int r1, const Eigen::MatrixXcd &op1);
  MPO OneBody_backward(const std::vector<int> & dims2, int r1, const Eigen::MatrixXcd &op1);
  MPO OneBody_forward_backward(const std::vector<int> & dims2, int r1, const Eigen::MatrixXcd &op1, const Eigen::MatrixXcd &op2);
  MPO OneBody_Sum(const std::vector<Eigen::MatrixXcd> & list);

  MPO TwoBody(const std::vector<int> & dims, int r1, int r2, 
                 const Eigen::MatrixXcd &op1, const Eigen::MatrixXcd &op2);
  MPO TwoBody_forward(const std::vector<int> & dims2, int r1, int r2, 
                 const Eigen::MatrixXcd &op1, const Eigen::MatrixXcd &op2);
  MPO TwoBody_backward(const std::vector<int> & dims2, int r1, int r2, 
                 const Eigen::MatrixXcd &op1, const Eigen::MatrixXcd &op2);
  MPO TwoBody_forward_backward(const std::vector<int> & dims2, int r1, int r2, 
                 const Eigen::MatrixXcd &op1, const Eigen::MatrixXcd &op2);

  MPO Product(const std::vector<Eigen::MatrixXcd> & list);

  void add(MPO & mpo, const MPO & other);
  MPO exp_Taylor(const MPO & other, const std::complex<double> & scale, const TruncatedSVD & trunc, int N_Taylor);
//TODO: forward & backward with single matrix (via SVD);

}//namespace MPO_op
}//namespace ACE
#endif 
