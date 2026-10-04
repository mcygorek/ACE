#ifndef ACE_SIMULATION_TEMPO_DEFINED_H
#define ACE_SIMULATION_TEMPO_DEFINED_H

#include "Parameters.hpp"
#include "ProcessTensorBuffer.hpp"
#include "FreePropagator.hpp"
#include "InitialState.hpp"
#include "OutputPrinter.hpp"
#include "TruncationLayout.hpp"
#include "ProcessTensorForwardList.hpp"

namespace ACE{

class Simulation_TEMPO{
public:
  bool print_timesteps;
  bool no_system_propagation;
  std::string print_dims_file;
  int print_dims_step;

  static void initialize_PTB(ProcessTensorBuffer &PTB, int n_mem, const DiagBB & diagBB, const Eigen::VectorXcd & rho_reduced);
  
  static std::vector<Eigen::MatrixXcd> initialize_b(int n_mem, DiagBB &diagBB, double dt);

  void single_step(ProcessTensorBuffer & PTB, Propagator &prop, const std::vector<Eigen::MatrixXcd> &b, const TimeGrid &tgrid, int n, const TruncationLayout & trunc_layout)const;

  Eigen::MatrixXcd run(Propagator &prop, DiagBB &diagBB,
           const Eigen::MatrixXcd & initial_rho, const TimeGrid &tgrid,
           OutputPrinter &printer, const TruncationLayout &trunc_layout)const;

  void setup(Parameters &param);

  Simulation_TEMPO(Parameters &param){
    setup(param);
  }
  Simulation_TEMPO(){
    Parameters param;
    setup(param);
  }

};


}//namespace
#endif
