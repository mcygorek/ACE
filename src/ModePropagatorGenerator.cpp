#include "ModePropagatorGenerator.hpp"
#include "ModePropagator.hpp"
#include "Potential1D.hpp"
#include "Parameters.hpp"
#include "Operators.hpp"
#include "MPG_Discretization.hpp"
#include "DummyException.hpp"

namespace ACE{

  std::string ModePropagatorGenerator::add_name(const std::string &str)const{
    return std::string(name()+"_"+str);
  }

  // get nr. of next mode; possibly skip
  int ModePropagatorGenerator::next(int k_prior)const{  
    int k=k_prior+1;
    if(k<0 || k>=get_N_modes()){
      return get_N_modes();
    }
    if(!skip_list[k])return k;
    return next(k);
  }

  std::vector<Eigen::MatrixXcd> ModePropagatorGenerator::get_env_ops(int k)const{
    std::vector<Eigen::MatrixXcd> mats;
    return mats;
  } 

  int ModePropagatorGenerator::get_mode_dim(int k)const{
    if(get_N_modes()<1){
      std::cerr<<"ModePropagatorGenerator::get_mode_dim: get_N_modes()<1! Not set up?"<<std::endl; 
      exit(1);
    }
    return get_bath_init(0).rows();
  }

  void ModePropagatorGenerator::setup_skip(Parameters &param){
    skip_list=std::vector<bool>(get_N_modes(), false);

    std::vector<size_t> skips=param.get_all_size_t(add_name("skip_mode"));
    for(size_t i=0; i<skips.size(); i++){
      int k=skips[i];
      if(k>=get_N_modes()){
        std::cerr<<add_name("skip_mode")<<" out of bounds: "<<k<<"/"<<get_N_modes()<<std::endl;
        exit(1);
      }
      skip_list[k]=true;
      skip_was_set=true;
    }
  }
  void ModePropagatorGenerator::zero_pad(int N_new){
    std::cerr<<"ModePropagatorGenerator::zero_pad: To be implemented in dependent classes!"<<std::endl;
    throw DummyException();
  }
  void ModePropagatorGenerator::setup_default(Parameters &param){
    set_N_modes(param.get_as_size_t(add_name("N_modes")));
    setup_skip(param);
  }
 

}//namespace

