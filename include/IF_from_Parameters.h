#ifndef IF_FROM_PARAMETERS_DEFINED_H
#define IF_FROM_PARAMETERS_DEFINED_H

#include "Parameters.h"
#include "InitialState.h"
#include "InfluenceFunctional_OD.h"
#include "InfluenceFunctional_Krylov.h"
#include "RankCompressor_Selector.h"


std::vector<Smart_Ptr<InfluenceFunctional_OD> > IF_from_Parameters(Parameters &param){

  std::vector<Smart_Ptr<InfluenceFunctional_OD> >IF;

  IF_TimeGrid tgrid(param);
  if(param.get_as_double("print_timegrid_info",true))tgrid.print_info();


  std::string read_PT=param.get_as_string("read_PT", "");
  std::string write_PT=param.get_as_string("write_PT", "");

  bool use_dict=param.get_as_bool("use_dict",false);
  double dict_zero=param.get_as_double("dict_zero", use_dict? 1e-12 : -1.);
  int factorization=param.get_as_int("factorization", 0);
  bool use_process_tensor=param.get_as_bool("use_process_tensor");


  int dim=2;
  if(param.is_specified("initial")){
    Eigen::MatrixXcd rho=InitialState(param);
    dim=rho.rows();
  }


  std::cout<<"Setting up IF"<<std::endl;
  RankCompressor_Ptr compressor=RankCompressor_Selector(param);


  if(read_PT!=""){
    if(use_process_tensor){
      std::cerr<<"'read_PT' cannot be used together with 'use_process_tensor'!"<<std::endl;
      exit(1);
    }
    std::cout<<"Reading process tensor file '"<<read_PT<<"'"<<std::endl;
    IF.push_back(new InfluenceFunctional_OD(read_PT)); 

  }else if(use_process_tensor){
    RealFunctionPtr SD=new RealFunction_Zero_Class();
    if(param.get_as_bool("use_bath",true))SD=new SpectralDensity_QD();

    double temperature=param.get_as_double("temperature",4);

    Operators2x2 op;
    Eigen::MatrixXcd couplings=param.get_as_operator("process_tensor_couplings",op.ketbra(1,1));
    DiagBB diagBB(couplings, SD, temperature);
     
    IF.push_back(new InfluenceFunctional_OD(tgrid, diagBB, compressor.ref()));

  }else{
    IF.push_back(new InfluenceFunctional_OD(tgrid, dim));
  }


  IF[0]->print_timesteps=param.get_as_bool("IF_print_timesteps",false);

  if(tgrid.use_rep){
    std::cout<<"rep_unit: "<<IF[0]->tgrid.rep_unit<<" n_rep: "<<IF[0]->tgrid.n_rep<<std::endl;


    IF[0]->compress_trafo_use_ortho=param.get_as_bool("compress_trafo_use_ortho",false);
  }

  if(param.get_as_size_t("Leads_N_modes",0)>0){
    std::cout<<"Calculating Influence functional for coupling to leads"<<std::endl;
    ModePropagatorGenerator *mpg=new ModePropagatorGenerator_Leads(param);
    IF[0]->add_modes(*mpg, *compressor, dict_zero);
    delete mpg;
  }

  if(param.get_as_size_t("RadiativeDecay_N_modes",0)>0){
    std::cout<<"Calculating Influence functional for radiative decay"<<std::endl;
    ModePropagatorGenerator *mpg=new ModePropagatorGenerator_RadiativeDecay(param);
    IF[0]->add_modes(*mpg, *compressor, dict_zero);
    delete mpg;
  }
  if(param.get_as_size_t("Superradiance_N_modes",0)>0){
    std::cout<<"Calculating Influence functional for coupling to photon modes"<<std::endl;
    ModePropagatorGenerator *mpg=new ModePropagatorGenerator_Superradiance(param);
    IF[0]->add_modes(*mpg, *compressor, dict_zero);
    delete mpg;
  }
  if(param.get_as_size_t("RandomSpin_N_modes",0)>0){
    std::cout<<"Calculating Influence functional for coupling random spin bath"<<std::endl;
    ModePropagatorGenerator *mpg=new ModePropagatorGenerator_RandomSpin(param);
    IF[0]->add_modes(*mpg, *compressor, dict_zero);
    delete mpg;
  }

  if(param.get_as_size_t("QDPhonon_N_modes",0)>0){
    std::cout<<"Calculating Influence functional for coupling to phonons"<<std::endl;
    ModePropagatorGenerator *mpg=new ModePropagatorGenerator_QDPhonon(param);
    IF[0]->add_modes(*mpg, *compressor, dict_zero);
    delete mpg;
  }
 

std::cout<<"Test dict: ";IF[0]->dict.print_beta(); std::cout<<std::endl;



  //Insert repeated units
  if(tgrid.use_rep){
    if(tgrid.rep_regularize){
      IF[0]->rep.regularize();
    }
#ifdef DEBUG_REP
    std::cout<<"insert: rep.M: "<<IF[0]->rep.M.dim_i<<" "<<IF[0]->rep.M.dim_d1<<" "<<IF[0]->rep.M.dim_d2<<std::endl;
    IF[0]->print_dims(std::cout);
#endif
    
    //do the actual insertion of the repeated units
    IF[0]->insert_rep();
  }



  if(write_PT!=""){
    IF[0]->write_binary(write_PT);
  }


  if(param.is_specified("multi_PT")){
    std::vector<std::string> sv=param.get_all_strings("multi_PT");
    for(size_t i=0; i<sv.size(); i++){
      std::cout<<"Reading process tensor file ["<<i+1<<"]: '"<<sv[i]<<"'"<<std::endl;
      IF.push_back(new InfluenceFunctional_OD(sv[i])); 
    }
  }

  
  if(param.get_as_size_t("system_cavity_QDPhonon_N_modes",0)>0){
    std::cout<<"Calculating Influence functional for coupling to phonons for system with cavity"<<std::endl;
    ModePropagatorGenerator *mpg=new ModePropagatorGenerator_system_cavity_QDPhonon(param);
    IF[0]->add_modes(*mpg, *compressor, dict_zero);
    delete mpg;
  }


  return IF;
}

#endif
