#ifndef ACE_GENERIC_SIMULATION_DEFINED_H
#define ACE_GENERIC_SIMULATION_DEFINED_H

#include "Parameters.hpp"
#include "OutputPrinter.hpp"
#include "Simulation_PT.hpp"
#include "ProcessTensorForwardList.hpp"
#include "InitialState.hpp"

namespace ACE{

class GenericSimulation{
public:
  Eigen::MatrixXcd initial_rho;
  int sysdim;
  FreePropagator fprop;
  OutputPrinter printer;
  Simulation_PT sim;
  int repeat_propagation;
  std::string write_final_densmat;
  
  inline int get_sysdim()const{return sysdim;}

  void setup(Parameters &param);

  Eigen::MatrixXcd run(const TimeGrid &tgrid, ProcessTensorForwardList &PT);

  GenericSimulation(Parameters &param){
    setup(param);
  }
};

}//namespace
#endif
