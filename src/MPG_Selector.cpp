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
    if(param.get_as_size_t("MultiSite_N_modes",0)>0
      || param.get_as_string("MultiSite_E_g_from_table")!=""){
      mpgs.push_back(std::make_shared<ModePropagatorGenerator_MultiSite>(param));
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
/*
    if(param.get_as_size_t("QDPhonon_N_modes",0)>0){
      mpgs.push_back(new ModePropagatorGenerator_QDPhonon(param));
    }
*/
    if(param.get_as_size_t("Lasing_N_modes",0)>0){
      mpgs.push_back(std::make_shared<ModePropagatorGenerator_Lasing>(param));
    }
    if(param.is_specified("add_single_mode")){
      std::vector<std::vector<std::string> > lines=param.get("add_single_mode");
      for(size_t r=0; r<lines.size(); r++){
        mpgs.push_back(std::make_shared<ModePropagatorGenerator_SingleMode>(lines[r]));
      }
    }
    if(param.is_specified("add_single_mode_from_file")){
      std::vector<std::vector<std::string> > lines=param.get("add_single_mode_from_file");
      for(size_t r=0; r<lines.size(); r++){
        mpgs.push_back(std::make_shared<ModePropagatorGenerator_SingleModeFromFile>(lines[r]));
      }
    }
    if(param.is_specified("interweave_modes")){
      std::vector<std::vector<std::string> > lines=param.get("interweave_modes");
      for(size_t r=0; r<lines.size(); r++){
        mpgs.push_back(std::make_shared<ModePropagatorGenerator_Interweave>(lines[r], param));
      }
    }
  }
  

}//namespace
