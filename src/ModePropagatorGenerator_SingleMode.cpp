#include "ModePropagatorGenerator_SingleMode.hpp"
#include "ModePropagatorGenerator.hpp"
#include "Parameters.hpp"
#include "Operators.hpp"
#include "Operators_Boson.hpp"
#include "ReadTable.hpp"
#include "Reader.hpp"

namespace ACE{

  void ModePropagatorGenerator_SingleMode::check_validity()const{
    int N=get_N();
    if(HE.rows()!=N*rho_init.rows()){
      std::cerr<<"SingleMode: Dimension of HE must be factorizable into N_sys*N_env, whereas rho_E_init must be of dimension N_env!"<<std::endl;
      std::cerr<<"Hint: Please check parameters 'add_single_mode [HE] [rho_E_init]'!"<<std::endl;
      exit(1);
    }
  }

  void ModePropagatorGenerator_SingleMode::setup(const Eigen::MatrixXcd & HE_,  const Eigen::MatrixXcd & init_){
    set_N_modes(1);
    HE=HE_;
    rho_init=init_;
    check_validity();
  }

  void ModePropagatorGenerator_SingleMode::setup(const std::vector<Eigen::MatrixXcd> & ops){
    set_N_modes(1);
    if(ops.size()<2){
      std::cerr<<"ModePropagatorGenerator_SingleMode: ops.size()<2!"<<std::endl;
      exit(1);
    }
    HE=ops[0];
    rho_init=ops[1];

std::cout<<"SingleMode: HE: "<<std::endl<<HE<<std::endl;

    envops.clear();
    for(size_t i=2; i<ops.size(); i++)envops.push_back(ops[i]);
    check_validity();
  }

  void ModePropagatorGenerator_SingleMode::setup(const std::vector<std::string> & str){
    std::vector<Eigen::MatrixXcd> ops;
    std::cout<<"Add single environment mode:"<<std::flush;
    for(size_t i=0; i<str.size(); i++){
      std::cout<<" '"<<str[i]<<"'"<<std::flush;
      ops.push_back(ReadExpression(str[i]));
    }
    std::cout<<std::endl;
    setup(ops);
  }

  ModePropagatorPtr ModePropagatorGenerator_SingleMode::getModePropagator(int k)const{
    if(k<0||k>=get_N_modes()){
      std::cerr<<"ModePropagatorGenerator_SingleMode: k<0||k>=get_N_modes()!"<<std::endl; 
      exit(1);
    }

    int sysdim=get_N();

    return ModePropagatorPtr(new ModePropagator(sysdim,get_bath_init(k),HE,get_env_ops(k)));
  }

}//namespace
