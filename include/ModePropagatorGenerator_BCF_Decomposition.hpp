#ifndef MODE_PROPAGATOR_GENERATOR_BCF_DECOMPOSITION_DEFINED_H
#define MODE_PROPAGATOR_GENERATOR_BCF_DECOMPOSITION_DEFINED_H

#include <vector>
#include <Eigen/Dense>
#include "ModePropagatorGenerator.hpp"
#include "BCF_Decomposition.hpp"

namespace ACE{

class ModePropagatorGenerator_BCF_Decomposition: public ModePropagatorGenerator{
public:

  std::vector<BCF_Decomposition::Term> terms; 
  
  int M;
  Eigen::MatrixXcd sysop;

  inline virtual std::string name()const{return std::string("BCF_Decomposition");}

  inline virtual int get_N()const{ 
    if(sysop.rows()<=2) return 2; 
    else return sysop.rows();
  }
  
  virtual EnvironmentOperators get_env_ops(int k)const;
  
  virtual void setup(Parameters &param);

  virtual ModePropagatorPtr getModePropagator(int k)const;

  virtual Eigen::MatrixXcd get_bath_init(int k)const;

  inline ModePropagatorGenerator_BCF_Decomposition(Parameters &param){
    setup(param);
  }
};

}//namespace
#endif
