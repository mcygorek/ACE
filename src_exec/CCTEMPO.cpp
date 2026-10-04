#include "Simulation_CCTEMPO.hpp"
#include "CCTEMPO.hpp"

using namespace ACE;

int main(int args, char** argv){

 try{
  Parameters param(args, argv, true);
  
  //int N_sites = CCTEMPO::get_N_sites(param);
  //std::vector<int> dims = CCTEMPO::get_dimensions(param);

  TimeGrid tgrid(param);
  TruncationLayout trunc_layout(param);

  // Set up MPO-form propagator
  MPO_Propagator mpo_prop(param);
  // Nodes of CTEMPO networks
  CCTEMPO_IF_Nodes nodes(param);
  //Initial state
  MPO initial=CCTEMPO::get_initial_state(param);

  //initialize comb-shaped state with teeth
  //Comb comb(param);

  Simulation_CCTEMPO sim(param);

  OutputPrinter_CCTEMPO outp(param);

  sim.run(mpo_prop, nodes, initial, tgrid, trunc_layout, outp);

 }catch (DummyException &e){
  return 1;
 }
#ifdef EIGEN_USE_MKL_ALL
  mkl_free_buffers();
#endif
  return 0;
}
 
