#ifndef ACE_SIMULATION_CCTEMPO_DEFINED_H
#define ACE_SIMULATION_CCTEMPO_DEFINED_H

#include "Parameters.hpp"
#include "FreePropagator.hpp"
#include "Simulation_CTEMPO.hpp"
#include "Comb.hpp"
#include "MPO_Propagator.hpp"
#include "MPO_ApplyOperators.hpp"
#include "TimeGrid.hpp"
#include "OutputPrinter_CCTEMPO.hpp"
#include "CCTEMPO_IF_Nodes.hpp"

namespace ACE{

class Simulation_CCTEMPO{
public:
  
  FreePropagator dummy_prop;
  Simulation_CTEMPO simC;
  bool use_symmetric_Trotter;
  int buffer_blocksize;
  std::string buffer_filename;

  Eigen::MatrixXcd SX_Hamiltonian_for_eigenstates;

  Eigen::MatrixXcd get_SX_Hamiltonian_for_eigenstates(double t);

  void run(MPO_Propagator &mpo_prop, CCTEMPO_IF_Nodes &nodes, const MPO & initial, const TimeGrid &tgrid, const TruncationLayout &trunc_layout, OutputPrinter_CCTEMPO &outp);

  void setup(Parameters &param);

  Simulation_CCTEMPO(Parameters &param){
    setup(param);
  }
  Simulation_CCTEMPO(){
    Parameters param;
    setup(param);
  }

};


}//namespace
#endif
