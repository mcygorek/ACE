#include "GenericSimulation.hpp"
#include "InitialState.hpp"
#include "BinaryReader.hpp"

namespace ACE{

void GenericSimulation::setup(Parameters &param){

  initial_rho=InitialState(param);
  sysdim=initial_rho.rows();

  fprop = FreePropagator(param, sysdim);

  printer.setup(param, sysdim);
  sim.setup(param);

  repeat_propagation = param.get_as_size_t("repeat_propagation",0);
  write_final_densmat = param.get_as_string("write_final_densmat","");
}

Eigen::MatrixXcd GenericSimulation::run(const TimeGrid &tgrid, ProcessTensorForwardList &PT){
  //run once
  Eigen::MatrixXcd state = sim.run(fprop, PT, initial_rho, tgrid, printer);

  //optionally: repeat
  TimeGrid tgrid2=tgrid; 
  for(size_t i=0; i<repeat_propagation; i++){
    std::cout<<"repeat: "<<i<<"/"<<repeat_propagation<<std::endl;
    tgrid2.ta=tgrid2.get_t_tot(); 
    state = sim.run(fprop, PT, state, tgrid2, printer);
  }

  //write final density matrix to file
  if(write_final_densmat!=""){
    binary_write_EigenMatrixXcd(write_final_densmat, state);
  }

  return state;
}

}//namespace
