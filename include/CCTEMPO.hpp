#pragma once
#ifndef CCTEMPO_DEFINED_H
#define CCTEMPO_DEFINED_H

#include "Parameters.hpp"
#include "FreePropagator.hpp"
#include "MPO.hpp"
#include "TruncatedSVD.hpp"

namespace ACE{
namespace CCTEMPO{

std::vector<int> get_dimensions(Parameters &param);
std::vector<int> get_dimensions_squared(Parameters &param);
int get_N_sites(Parameters &param);

TruncatedSVD get_SysProp_trunc(Parameters &param);
TruncatedSVD get_SysPropMPO_trunc(Parameters &param);

inline bool get_symmetric_Trotter(Parameters &param){ 
  return param.get_as_bool("use_symmetric_Trotter", true);
}

FreePropagator get_SX_FreePropagator(Parameters &param);
Eigen::MatrixXcd get_SX_Hamiltonian(Parameters &param, double t=0);
Eigen::MatrixXcd get_SX_Hamiltonian_from_individual_terms(Parameters &param);

MPO get_initial_state(Parameters &param);

void add_HilbertSpace_Propagator_CCTEMPO(std::vector<MPO> & sys_prop, Parameters &param, double t=0);

void add_collective_decay(std::vector<MPO> & sys_prop, Parameters & param);
std::vector<MPO> get_sys_prop_mpo(Parameters &param, double t=0);

}//namespace CCTEMPO
}//namespace ACE
#endif
