#ifndef MODE_PROPAGATOR_GENERATOR_INTERWEAVE_DEFINED_H
#define MODE_PROPAGATOR_GENERATOR_INTERWEAVE_DEFINED_H

#include "ModePropagatorGenerator.hpp"
#include <iostream>
#include "DummyException.hpp"

namespace ACE{

class ModePropagatorGenerator_Interweave: public ModePropagatorGenerator{
public:
  std::vector<std::shared_ptr<ModePropagatorGenerator> > mpgs;
  //translation: joint mode index -> (which generator, mode index within generator)
  std::vector<std::pair<int,int> > which;
 

  //Name to identify sets of parameters for a give MPG:
  virtual std::string name()const;
  
  //'label' for mode index [used for printing], e.g., rescaling to frequency
  //inline virtual double k_label(int k)const{ return k; }

  virtual int get_N()const; // system dimension 

  inline virtual int get_N_modes()const{  // get number of modes
    return which.size();
  }
  inline virtual void set_N_modes(int Nmodes){  // set number of modes
    std::cerr<<"Cannot 'set_N_modes' in ModePropagatorGenerator_Interweave!"<<std::endl;
    throw DummyException();
    //skip_list.resize(Nmodes, false);
  }

  virtual std::vector<Eigen::MatrixXcd> get_env_ops(int k)const;
   
  virtual Eigen::MatrixXcd get_bath_init(int k)const;

  virtual int get_mode_dim(int k=0)const;


  void check_exists(int k, const std::string &context)const;
  void print_which(std::ostream &os=std::cout)const;

  void setup(const std::vector<std::string> & str, Parameters &param);


  virtual ModePropagatorPtr getModePropagator(int k)const;

  inline ModePropagatorGenerator_Interweave(const std::vector<std::string> & str, Parameters &param){
    setup(str, param);
  }
  inline ModePropagatorGenerator_Interweave(){
  }
  virtual ~ModePropagatorGenerator_Interweave(){}

};

}//namespace

#endif
