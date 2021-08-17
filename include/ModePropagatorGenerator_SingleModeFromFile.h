#ifndef MODE_PROPAGATOR_GENERATOR_SINGLEMODE_FROM_FILE_DEFINED_H
#define MODE_PROPAGATOR_GENERATOR_SINGLEMODE_FROM_FILE_DEFINED_H

#include "ModePropagatorGenerator.h"
#include "Parameters.h"
#include "Operators.h"
#include "Operators_Boson.h"
#include "ReadTable.h"

class ModePropagatorGenerator_SingleModeFromFile: public ModePropagatorGenerator{
public:

  ModePropagatorPtr mpp; 
  Eigen::MatrixXcd rho_init;

  std::vector<Eigen::MatrixXcd> envops;

  virtual std::string name()const{return std::string("SingleModeFromFile");}

  virtual int get_N()const{ return mpp->get_N_system(); }

  void setup(const std::string &file, const std::vector<Eigen::MatrixXcd> & ops){
    set_N_modes(1);
    if(ops.size()<1){
      std::cerr<<"ModePropagatorGenerator_SingleModeFromFile: ops.size()<1!"<<std::endl;
      exit(1);
    }

    rho_init=ops[0];
    envops.clear();
    for(size_t i=1; i<ops.size(); i++)envops.push_back(ops[i]);

    Parameters param2;
    param2.add_from_file(file);
    mpp=new ModePropagator(FreePropagator(param2), rho_init, envops);

  }
  void setup(const std::vector<std::string> & str){
    if(str.size()<2){
      std::cerr<<"ModePropagatorGenerator_SingleModeFromFile: needs at least 2 arguments: filename bath_init!"<<std::endl;
      exit(1);
    }
    std::vector<Eigen::MatrixXcd> ops;
//    std::cout<<"Add single environment mode:"<<std::flush;
    for(size_t i=1; i<str.size(); i++){
      std::cout<<" '"<<str[i]<<"'"<<std::flush;
      ops.push_back(ReadExpression(str[i]));
    }
//    std::cout<<std::endl;
    setup(str[0],ops);
  }

  virtual std::vector<Eigen::MatrixXcd> get_env_ops() const{
    return envops;
  }
  virtual Eigen::MatrixXcd get_bath_init(int k)const{
    return rho_init;
  }
  virtual ModePropagatorPtr getModePropagator(int k)const{
    if(k<0||k>=get_N_modes()){
      std::cerr<<"ModePropagatorGenerator_SingleMode: k<0||k>=get_N_modes()!"<<std::endl; 
      exit(1);
    }

    return mpp;
  }

  ModePropagatorGenerator_SingleModeFromFile(const std::string &file, const std::vector<Eigen::MatrixXcd> &ops){
    setup(file, ops);
  }
  ModePropagatorGenerator_SingleModeFromFile(const std::vector<std::string> &svec){
    setup(svec);
  }
};

#endif
