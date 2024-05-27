#include "ModePropagatorGenerator_BEC.hpp"
#include "ModePropagatorGenerator.hpp"
#include "Parameters.hpp"
#include "Operators.hpp"
#include "otimes.hpp"
#include "Operators_Boson.hpp"
#include "ReadTemperature.hpp"

namespace ACE{

void ModePropagatorGenerator_BEC::setup(Parameters &param){
  E_g.setup(param, name());
  set_N_modes(E_g.N);
  setup_skip(param);

  M=param.get_as_size_t(add_name("M"), 2);
  N=param.get_as_size_t(add_name("N"), 2);
  temperature=readTemperature(param,name());

} 

ModePropagatorPtr ModePropagatorGenerator_BEC::getModePropagator(int k)const{
  if(k<0||k>=get_N_modes()){
    std::cerr<<"ModePropagatorGenerator_BEC: k<0||k>=get_N_modes()!"<<std::endl; 
    exit(1);
  }
  
  Eigen::MatrixXcd Htot=Eigen::MatrixXcd::Zero(N*M,N*M);  

  return ModePropagatorPtr(new ModePropagator(get_N(),get_bath_init(k),Htot,get_env_ops(k)));
}

Eigen::MatrixXcd ModePropagatorGenerator_BEC::get_bath_init(int k)const{
  return Boson_Equilibrium(Operators_Boson::n(M)*get_E(k), temperature);
}

}//namespace
