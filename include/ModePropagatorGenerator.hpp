#ifndef MODE_PROPAGATOR_GENERATOR_DEFINED_H
#define MODE_PROPAGATOR_GENERATOR_DEFINED_H

#include <vector>
#include "Eigen_fwd.hpp"
#include "ModePropagator.hpp"
#include "Potential1D.hpp"
#include "Parameters.hpp"
#include "Operators.hpp"
#include "EnvironmentOperators.hpp"
#include "MPG_Discretization.hpp"

namespace ACE{

class ModePropagatorGenerator{
public:
  std::vector<bool> skip_list; //Indicates if the calculation of a certain mode is to be skipped
  bool skip_was_set;  // terminate with an error if set_N_modes is called after skipped modes were set. NOTE: Currently not used!

  
  //Name to identify sets of parameters for a give MPG:
  virtual std::string name()const=0;
  virtual std::string add_name(const std::string &str)const;
  

  //'label' for mode index [used for printing], e.g., rescaling to frequency
  inline virtual double k_label(int k)const{ return k; }

  virtual int get_N()const=0; // system dimension 

  inline virtual int get_N_modes()const{  // get number of modes
    return skip_list.size();
  }
  inline virtual void set_N_modes(int Nmodes){  // set number of modes
    skip_list.resize(Nmodes, false);
  }

  virtual int next(int k_prior)const;  // get nr. of next mode; possibly skip
  
  inline virtual int first()const{return next(-1);}

  inline virtual bool was_set_up()const{return get_N_modes()>0; }

  virtual EnvironmentOperators get_env_ops(int k) const;
   
  virtual Eigen::MatrixXcd get_bath_init(int k)const=0;

  virtual int get_mode_dim(int k=0)const;

  virtual void zero_pad(int N_new);

  void setup_skip(Parameters &param);
  
  void setup_default(Parameters &param);

  virtual ModePropagatorPtr getModePropagator(int k)const=0;

  inline ModePropagatorGenerator(){
    skip_was_set=false;
  }
  virtual ~ModePropagatorGenerator(){}

};

}//namespace

#endif
