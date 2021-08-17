#ifndef MODE_PROPAGATOR_GENERATOR_DEFINED_H
#define MODE_PROPAGATOR_GENERATOR_DEFINED_H

#include "ModePropagator.h"
#include "Parameters.h"
#include "Operators.h"
#include "MPG_Discretization.h"

class ModePropagatorGenerator{
public:
  std::vector<bool> skip_list; //Indicates if the calculation of a certain mode is to be skipped
  bool skip_was_set;  // terminate with an error if set_N_modes is called after skipped modes were set. NOTE: Currently not used!

  
  //Name to identify sets of parameters for a give MPG:
  virtual std::string name()const=0;
  virtual std::string add_name(const std::string &str)const{
    return std::string(name()+"_"+str);
  }

  //'label' for mode index [used for printing], e.g., rescaling to frequency
  virtual double k_label(int k)const{ return k; }

  virtual int get_N()const=0; // system dimension 

  virtual int get_N_modes()const{  // get number of modes
    return skip_list.size();
  }
  virtual void set_N_modes(int Nmodes){  // set number of modes
/*
    if(skip_was_set){
      std::cerr<<"ModePropagatorGenerator: set_N_modes was called after skipping of certain modes was set!"<<std::endl;
      exit(1);
    }
    skip_list.clear();
*/
    skip_list.resize(Nmodes, false);
  }

  virtual int next(int k_prior)const{  // get nr. of next mode; possibly skip
    int k=k_prior+1;
    if(k<0 || k>=get_N_modes()){
      return get_N_modes();
    }
    if(!skip_list[k])return k;
    return next(k);
  }
  virtual int first()const{return next(-1);}

  virtual bool was_set_up()const{return get_N_modes()>0; }

  virtual std::vector<Eigen::MatrixXcd> get_env_ops() const{
    std::vector<Eigen::MatrixXcd> mats;
    return mats;
  } 
  virtual Eigen::MatrixXcd get_bath_init(int k)const=0;

  virtual int get_mode_dim(int k=0)const{
    if(get_N_modes()<1){
      std::cerr<<"ModePropagatorGenerator::get_mode_dim: get_N_modes()<1! Not set up?"<<std::endl; 
      exit(1);
    }
    return get_bath_init(0).rows();
  }

  void setup_skip(Parameters &param){
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
  void setup_default(Parameters &param){
    set_N_modes(param.get_as_size_t(add_name("N_modes")));
    setup_skip(param);
  }
 

  virtual ModePropagatorPtr getModePropagator(int k)const=0;

  ModePropagatorGenerator(){
    skip_was_set=false;
  }
  virtual ~ModePropagatorGenerator(){}

};




#include "ModePropagatorGenerator_SingleMode.h"
#include "ModePropagatorGenerator_SingleModeFromFile.h"
#include "ModePropagatorGenerator_Fermion.h"
#include "ModePropagatorGenerator_Boson.h"
#include "ModePropagatorGenerator_Potential1D.h"
#include "ModePropagatorGenerator_MultiSite.h"
#include "ModePropagatorGenerator_QDPhonon.h"
#include "ModePropagatorGenerator_RandomSpin.h"

#endif
