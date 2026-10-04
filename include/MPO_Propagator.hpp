#pragma once
#ifndef ACE_MPO_PROPAGATOR_DEFINED_H
#define ACE_MPO_PROPAGATOR_DEFINED_H

#include "MPO_op.hpp"
#include "FreePropagator.hpp"
#include "Parameters.hpp"
#include "TruncatedSVD.hpp"
#include "TimeGrid.hpp"
#include "Comb.hpp"
#include "MPO_ApplyOperators.hpp"


namespace ACE{

/** Logic: If there are multi-site Hamiltonians, include the time-independent
           single-site terms into the multi-site Hamiltonian.
           The multi-site Hamiltonians for forward and backward propagation
           are treated separately, and are exponentiated using a given Taylor
           order.
           If no multi-site Hamiltonian is given or if the single-site 
           Hamiltonian has a time-dependent contribution, the contributions
           that are not included in the multi-site term are exponentiated 
           separately using the corresponding fprop.
           This corresponds to an additional Trotter splitting. 
           Also collective Lindblad terms are exponentiated separately
           and Trotter-split from the Hamiltonian terms.
*/

class MPO_Propagator{
public:
  std::vector<std::shared_ptr<FreePropagator> > fprop;
  std::vector<MPO> sys_prop_mpo;
  MPO_ApplyOperators apply_ops;

  bool use_symmetric_Trotter;

  
  void propagate_comb(Comb &comb, const TimeGrid &tgrid, int step, const TruncatedSVD &trunc);
  //for second part of symmetric Trotter splitting:
  void propagate_comb_second(Comb &comb, const TimeGrid &tgrid, int step, const TruncatedSVD &trunc);

  void setup(Parameters &param);

  MPO_Propagator(Parameters &param){setup(param);}
  MPO_Propagator(){}
  ~MPO_Propagator(){}
};
}
#endif
