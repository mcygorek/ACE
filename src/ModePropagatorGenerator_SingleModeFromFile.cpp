#include "ModePropagatorGenerator_SingleModeFromFile.hpp"
#include "ModePropagatorGenerator.hpp"
#include "Parameters.hpp"
#include "Operators.hpp"
#include "Operators_Boson.hpp"
#include "Reader.hpp"

namespace ACE{


  void ModePropagatorGenerator_SingleModeFromFile::setup(Parameters &param2, const std::vector<Eigen::MatrixXcd> & ops){
    set_N_modes(1);
    if(ops.size()<1){
      std::cerr<<"ModePropagatorGenerator_SingleModeFromFile: ops.size()<1!"<<std::endl;
      exit(1);
    }

    envops.clear();
    for(size_t i=1; i<ops.size(); i++)envops.push_back(ops[i]);

    mpp=std::make_shared<ModePropagator>(FreePropagator(param2), ops[0], envops);
  }
  void ModePropagatorGenerator_SingleModeFromFile::setup(const std::string &file, const std::vector<Eigen::MatrixXcd> & ops){
    Parameters param2;
    param2.add_from_file(file);
    setup(param2, ops);
  }
  

  void ModePropagatorGenerator_SingleModeFromFile::setup(const std::vector<std::string> & str){
    if(str.size()<2){
      std::cerr<<"ModePropagatorGenerator_SingleModeFromFile: needs at least 2 arguments: filename bath_init!"<<std::endl;
      exit(1);
    }
    std::vector<Eigen::MatrixXcd> ops;
    for(size_t i=1; i<str.size(); i++){
      std::cout<<" '"<<str[i]<<"'"<<std::flush;
      ops.push_back(ReadExpression(str[i]));
    }
    setup(str[0],ops);
  }

  ModePropagatorPtr ModePropagatorGenerator_SingleModeFromFile::get_ModePropagator(int k)const{
    if(k<0||k>=get_N_modes()){
      std::cerr<<"ModePropagatorGenerator_SingleMode: k<0||k>=get_N_modes()!"<<std::endl; 
      exit(1);
    }

    return mpp;
  }

}//namespace
