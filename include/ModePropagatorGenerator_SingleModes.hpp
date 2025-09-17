#pragma once
#ifndef MODE_PROPAGATOR_GENERATOR_SINGLEMODES_DEFINED_H
#define MODE_PROPAGATOR_GENERATOR_SINGLEMODES_DEFINED_H

#include "ModePropagatorGenerator.hpp"
#include "ModePropagatorGenerator_SingleModeFromFile.hpp"
#include <vector>
#include <Eigen/Dense>

namespace ACE{


class ModePropagatorGenerator_SingleModes: public ModePropagatorGenerator{
public:

  std::vector<ModePropagatorPtr> modes;

  inline virtual std::string name()const{return std::string("SingleModes");}

  inline virtual int get_N()const{ 
    if(modes.size()<1)return 0;
    else return modes[0]->get_N_system(); 
  }
  virtual void zero_pad(int N_new);

  void setup(Parameters &param);
  void check_validity()const;
  void add_single_mode(const std::vector<Eigen::MatrixXcd> & ops);
  void add_single_mode(const std::vector<std::string> & str);
  void add_single_mode_from_file(const std::string &file, const std::vector<Eigen::MatrixXcd> & ops);

  void add_single_mode_from_file(const std::vector<std::string> & str);

  void add_single_mode(std::vector<std::shared_ptr<ModePropagator> > & list);

  inline void add_single_mode(std::shared_ptr<ModePropagator> & SM){
    std::vector<std::shared_ptr<ModePropagator> > list(1, SM);
    add_single_mode(list);
  }

  virtual EnvironmentOperators get_env_ops(int k) const;
  virtual Eigen::MatrixXcd get_bath_init(int k)const;
  virtual ModePropagatorPtr getModePropagator(int k)const;

  ModePropagatorGenerator_SingleModes(Parameters &param){
    setup(param);
  }
  inline ModePropagatorGenerator_SingleModes(std::vector<std::shared_ptr<ModePropagator> > & list){
    add_single_mode(list);
  }
};

}//namespace
#endif
