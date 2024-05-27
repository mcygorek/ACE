#include "ModePropagatorGenerator_Lasing.hpp"
#include "Operators.hpp"
#include "Operators_Boson.hpp"
#include "otimes.hpp"
#include "InitialState.hpp"

namespace ACE{

std::vector<Eigen::MatrixXcd> ModePropagatorGenerator_Lasing::get_env_ops(int k)const{
  std::vector<Eigen::MatrixXcd> mats(3);
  mats[0]=0.5*sigma_x();
  mats[1]=0.5*sigma_y();
  mats[2]=0.5*sigma_z();
  return mats;
}

void ModePropagatorGenerator_Lasing::setup(Parameters &param){
  E_g.setup(param, name());
  set_N_modes(param.get_as_size_t(add_name("N_modes")));
  setup_skip(param);

  N_system = param.get_as_int(add_name("N_system"), -1);
  if(N_system<2){
    Eigen::MatrixXcd rho_init = InitialState(param);
    N_system = rho_init.rows();
  } 
  g = param.get_as_double_check(add_name("g"));
  g_prime = param.get_as_double(add_name("g_prime"),0);
  Gamma_up = param.get_as_double(add_name("Gamma_up"), 0);
  Gamma_down = param.get_as_double(add_name("Gamma_down"), 0);

  std::cout<<"Parameters for ModePropagatorGenerator_Lasing:";
  std::cout<<" N_system="<<N_system<<" g="<<g<<" g'="<<g_prime;
  std::cout<<" Gamma_up="<<Gamma_up<<" Gamma_down="<<Gamma_down<<std::endl;
}

ModePropagatorPtr ModePropagatorGenerator_Lasing::getModePropagator(int k)const{
  if(k<0||k>=get_N_modes()){
    std::cerr<<"ModePropagatorGenerator_Lasing: k<0||k>=get_N_modes()!"<<std::endl; 
    exit(1);
  }
  Operators op2(2);
  
  Eigen::MatrixXcd H = hbar_in_meV_ps * E_g.get_omega(k) * 
         otimes(Operators_Boson::id(N_system), op2.ketbra(1,1));

  H += hbar_in_meV_ps * E_g.get_g(k) * (
         otimes(Operators_Boson::a(N_system),sigma_plus()) 
       + otimes(Operators_Boson::adagger(N_system), sigma_minus()));

  H += hbar_in_meV_ps * g_prime * (
         otimes(Operators_Boson::a(N_system),sigma_minus()) 
       + otimes(Operators_Boson::adagger(N_system), sigma_plus()));

  FreePropagator fprop;
  fprop.add_Hamiltonian(H);
  fprop.add_Lindblad(Gamma_up, otimes(Operators_Boson::id(N_system),sigma_plus()));
  fprop.add_Lindblad(Gamma_down, otimes(Operators_Boson::id(N_system),sigma_minus()));
  
  return ModePropagatorPtr(new ModePropagator(fprop, get_bath_init(k), get_env_ops(k)));
}

Eigen::MatrixXcd ModePropagatorGenerator_Lasing::get_bath_init(int k)const{
  return Operators(2).ketbra(0,0);
}


}//namespace
