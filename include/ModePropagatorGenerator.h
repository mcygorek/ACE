#ifndef MODE_PROPAGATOR_GENERATOR_DEFINED_H
#define MODE_PROPAGATOR_GENERATOR_DEFINED_H

#include "ModePropagator.h"
#include "Parameters.h"
#include "Operators.h"

class ModePropagatorGenerator{
public:

  virtual int get_N()const=0;
  virtual int get_N_modes()const =0;
  virtual bool was_set_up()const{return get_N_modes()>0; }
  virtual std::vector<Eigen::MatrixXcd> get_env_ops() const{
    std::vector<Eigen::MatrixXcd> mats;
    return mats;
  } 
  virtual Eigen::MatrixXcd get_bath_init(int k)const=0;
  virtual ModePropagatorPtr getModePropagator(int k)const=0;
  virtual ~ModePropagatorGenerator(){}
};

#include "ModePropagatorGenerator_Leads.h"
#include "ModePropagatorGenerator_RadiativeDecay.h"
#include "ModePropagatorGenerator_QDPhonon.h"
#include "ModePropagatorGenerator_system_cavity_QDPhonon.h"
#include "ModePropagatorGenerator_Superradiance.h"
//#include "ModePropagatorGenerator_Superradiance3d.h"
#include "ModePropagatorGenerator_RandomSpin.h"
#endif
