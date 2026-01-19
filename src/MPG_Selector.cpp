#include "MPG_Selector.hpp"
#include "Parameters.hpp"
#include "ModePropagatorGeneratorList.hpp"

namespace ACE{

  void MPG_Selector::setup(Parameters &param){
    bool use_Gaussian=param.get_as_bool("use_Gaussian");
    std::string Gaussian_prefix=param.get_as_string("Gaussian_prefix","Boson");

    if(!(use_Gaussian && Gaussian_prefix=="Boson")  && 
       (param.get_as_size_t("Boson_N_modes",0)>0
        || param.get_as_string("Boson_E_g_from_table")!="") ){
      mpgs.push_back(std::make_shared<ModePropagatorGenerator_Boson>(param));
    }
    if(param.get_as_size_t("Potential1D_N_modes",0)>0
      || param.get_as_string("Potential1D_E_g_from_table")!=""){
      mpgs.push_back(std::make_shared<ModePropagatorGenerator_Potential1D>(param));
    }
    if(param.get_as_size_t("Fermion_N_modes",0)>0
      || param.get_as_string("Fermion_E_g_from_table")!=""){
      mpgs.push_back(std::make_shared<ModePropagatorGenerator_Fermion>(param));
    }
    if(param.get_as_size_t("RandomSpin_N_modes",0)>0){
      mpgs.push_back(std::make_shared<ModePropagatorGenerator_RandomSpin>(param));
    }
    if(param.is_specified("add_single_mode")||param.is_specified("add_single_mode_from_file")){
      mpgs.push_back(std::make_shared<ModePropagatorGenerator_SingleModes>(param));
    }
  }
  

}//namespace
