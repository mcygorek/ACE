#ifndef MODE_PROPAGATOR_GENERATOR_SINGLEMODE_FROM_FILE_DEFINED_H
#define MODE_PROPAGATOR_GENERATOR_SINGLEMODE_FROM_FILE_DEFINED_H

#include "ModePropagatorGenerator.hpp"
#include <vector>
#include <Eigen/Dense>

namespace ACE{

class ModePropagatorGenerator_SingleModeFromFile: public ModePropagatorGenerator{
public:

  ModePropagatorPtr mpp; 
  Eigen::MatrixXcd rho_init;

  std::vector<Eigen::MatrixXcd> envops;

  inline virtual std::string name()const{return std::string("SingleModeFromFile");}

  inline virtual int get_N()const{ return mpp->get_N_system(); }

  void setup(const std::string &file, const std::vector<Eigen::MatrixXcd> & ops);

  void setup(const std::vector<std::string> & str);

  inline virtual std::vector<Eigen::MatrixXcd> get_env_ops(int k) const{
    return envops;
  }
  inline virtual Eigen::MatrixXcd get_bath_init(int k)const{
    return rho_init;
  }
  virtual ModePropagatorPtr getModePropagator(int k)const;

  ModePropagatorGenerator_SingleModeFromFile(const std::string &file, const std::vector<Eigen::MatrixXcd> &ops){
    setup(file, ops);
  }
  ModePropagatorGenerator_SingleModeFromFile(const std::vector<std::string> &svec){
    setup(svec);
  }
};

}//namespace
#endif
