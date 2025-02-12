#include "PCH.hpp"
#include "ProcessTensorForwardList.hpp"
#include "DummyException.hpp"
#include "Parameters.hpp"
#include "EDM_Simulation.hpp"
#include "LiouvilleTools.hpp"
#include <Eigen/Core>

using namespace ACE;

int main(int args, char** argv){

  Parameters param(args,argv,true);

  EDM_Simulation sim(param);
  EDM_State state=EDM_Initial_State(param);

  std::cout<<"Initial state:"<<std::endl;
  state.print();
  std::cout<<"Closures:"<<std::endl;
  state.print_closures();

  for(int r=0; r<state.Ldim.size(); r++){
    std::cout<<std::endl;
    std::cout<<"Trace over "<<r<<":"<<std::endl;
    EDM_State tmp=state.trace_over(r); 
    tmp.print();
  }

  for(int r=0; r<state.Ldim.size(); r++){
    std::cout<<std::endl;
    std::cout<<"get_reduced("<<r<<"):"<<std::endl;
    std::cout<<L_Vector_to_H_Matrix(state.get_reduced(r))<<std::endl;
  }

  return 0;
}

