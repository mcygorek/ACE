#include "PCH.hpp"
#include "InitialState.hpp"
#include "BinaryReader.hpp"

namespace ACE{

  void InitialState::setup(Parameters &param){
    initialize();
    bool print_initial=param.get_as_bool("print_initial",false);
    std::string str=param.get_as_string("initial","");
    std::string read_densmat=param.get_as_string("initial_read_densmat","");
   
    if(read_densmat!=""){
      rho=binary_read_EigenMatrixXcd(read_densmat);
    }else if(str==""||str=="TLS_GS"){
      rho=Operators(2).ketbra(0,0);
      if(print_initial)std::cout<<"Using initial state: 2LS ground state"<<std::endl;
    }else if(str=="TLS_EX"||str=="EX"){
      rho=Operators(2).ketbra(1,1);
      if(print_initial)std::cout<<"Using initial state: 2LS excited"<<std::endl;
    }else if(str=="TLS_COH"){
      rho=Eigen::MatrixXcd(2,2);
      rho<< 0.5, 0.5, 0.5, 0.5;
      if(print_initial)std::cout<<"Using initial state: 2LS coherent"<<std::endl;
    }else{
      if(print_initial)std::cout<<"Using explicit initial state: '"<<str<<"'"<<std::endl;
      rho=param.get_as_operator("initial");
      if(print_initial)std::cout<<rho<<std::endl;
    }
  }


}//namespace
