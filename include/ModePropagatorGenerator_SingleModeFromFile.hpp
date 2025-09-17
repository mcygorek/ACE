#ifndef MODE_PROPAGATOR_GENERATOR_SINGLEMODE_FROM_FILE_DEFINED_H
#define MODE_PROPAGATOR_GENERATOR_SINGLEMODE_FROM_FILE_DEFINED_H

#include "ModePropagatorGenerator.hpp"
#include <vector>
#include <Eigen/Dense>
#include "DummyException.hpp"

namespace ACE{

class ModePropagatorGenerator_SingleModeFromFile: public ModePropagatorGenerator{
public:

  ModePropagatorPtr mpp; 

  std::vector<Eigen::MatrixXcd> envops;

  inline virtual std::string name()const{return std::string("SingleModeFromFile");}

  inline virtual int get_N()const{ 
    if(!mpp){
      std::cerr<<"SingleModeFromFile::get_N: mpp not set!"<<std::endl;
      throw DummyException();
    }
    return mpp->get_N_system(); 
  }

  void setup(Parameters &param2, const std::vector<Eigen::MatrixXcd> & ops);
  void setup(const std::string &file, const std::vector<Eigen::MatrixXcd> & ops);

  void setup(const std::vector<std::string> & str);

  inline virtual EnvironmentOperators get_env_ops(int k) const{
    return EnvironmentOperators(envops);
  }
  inline virtual Eigen::MatrixXcd get_bath_init(int k)const{
    if(!mpp){
      std::cerr<<"SingleModeFromFile::get_bath_init: mpp not set!"<<std::endl;
      throw DummyException();
    }
    return mpp->bath_init;
  }
  virtual ModePropagatorPtr getModePropagator(int k)const;

  ModePropagatorGenerator_SingleModeFromFile(const std::string &file, const std::vector<Eigen::MatrixXcd> &ops){
    setup(file, ops);
  }
  ModePropagatorGenerator_SingleModeFromFile(const std::string &file, const Eigen::MatrixXcd &initial){
    std::vector<Eigen::MatrixXcd> ops(1, initial);
    setup(file, ops);
  }
  ModePropagatorGenerator_SingleModeFromFile(Parameters &param2, const Eigen::MatrixXcd &initial){
    std::vector<Eigen::MatrixXcd> ops(1, initial);
    setup(param2, ops);
  }
  ModePropagatorGenerator_SingleModeFromFile(const std::vector<std::string> &svec){
    setup(svec);
  }
};

}//namespace
#endif
