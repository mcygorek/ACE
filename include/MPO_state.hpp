#ifndef MPO_STATE_DEFINED_H_
#define MPO_STATE_DEFINED_H_

#include "MPO.hpp"

namespace ACE{
namespace MPO_state{  //features with clear interpretation of MPO as operator

  std::vector<Eigen::VectorXcd> get_cap_d2(const MPO & mpo);
  std::vector<Eigen::VectorXcd> get_cap_d4(const MPO & mpo);
  std::complex<double> evaluate(const MPO & mpo, const std::vector<int> &list);

  MPO from_SingleExcitationManifold(Eigen::MatrixXcd &rho, std::vector<int> dims=std::vector<int>());

  std::complex<double> dot(const MPO & mpo1, const MPO & mpo2);
} //namespace MPO_state
} //namespace ACE
#endif
