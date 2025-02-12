#ifndef ACE_SIMULATION_PT_DEFINED_H
#define ACE_SIMULATION_PT_DEFINED_H

#include "Parameters.hpp"
#include "TimeGrid.hpp"
#include "FreePropagator.hpp"
#include "ProcessTensorForwardList.hpp"
#include "OutputPrinter.hpp"
#include "InitialState.hpp"

namespace ACE{

class Simulation_PT{
public:
  bool print_timesteps;
  bool print_final_maxdim;
  bool propagate_alternate;
  bool use_symmetric_Trotter;

  //For Transfer Tensors:
  bool use_TT;
  int TT_n_from;
  int TT_n_mem;
  bool use_LT;
  std::string TT_print_norms;
  std::string ME_print_rates;
  std::string ME_print_L;
  

  static void propagate_system(Eigen::MatrixXcd & state, Propagator &prop, double t, double dt);


  void propagate_state(Eigen::MatrixXcd &state, int n, const TimeGrid &tgrid, Propagator &prop, ProcessTensorForwardList &PT)const;

  Eigen::MatrixXcd run_std(Propagator &prop, ProcessTensorForwardList &PT,
           const Eigen::MatrixXcd & initial_rho, const TimeGrid &tgrid,
           OutputPrinter &printer);

//  Eigen::MatrixXcd run_TT(Propagator &prop, ProcessTensorForwardList &PT,
//           const Eigen::MatrixXcd & initial_rho, const TimeGrid &tgrid,
//           OutputPrinter &printer);

  Eigen::MatrixXcd run(Propagator &prop, ProcessTensorForwardList &PT,
           const Eigen::MatrixXcd & initial_rho, const TimeGrid &tgrid,
           OutputPrinter &printer);

  inline void run_(  //used for pybind11
           FreePropagator &prop, ProcessTensorForwardList &PT,
           const InitialState & initial, const TimeGrid &tgrid,
           OutputPrinter &printer){
    run(prop, PT, initial.rho, tgrid, printer); 
  }

  void setup(Parameters &param);

  Simulation_PT(Parameters &param){
    setup(param);
  }
  Simulation_PT(){ 
    Parameters param;
    setup(param);
  }
};
}//namespace
#endif
