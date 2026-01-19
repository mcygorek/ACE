#include "ModePropagatorGenerator_SingleModes.hpp"
#include "DummyException.hpp"
#include "ReadExpression.hpp"

namespace ACE{


void ModePropagatorGenerator_SingleModes::setup(Parameters &param){
  if(param.is_specified("add_single_mode")){
    std::vector<std::vector<std::string> > lines=param.get("add_single_mode");
    for(size_t r=0; r<lines.size(); r++){
      add_single_mode(lines[r]);
    }
  }
  if(param.is_specified("add_single_mode_from_file")){
    std::vector<std::vector<std::string> > lines=param.get("add_single_mode_from_file");
    for(size_t r=0; r<lines.size(); r++){
      add_single_mode_from_file(lines[r]);
    }
  }

  std::cout<<"SingleModes: modes.size()="<<modes.size()<<" skip_list.size()="<<skip_list.size()<<std::endl;
}
void ModePropagatorGenerator_SingleModes::zero_pad(int N_new){ 
  int NS=get_N();
  ModePropagatorPtr zeromode(new ModePropagator(NS, Eigen::MatrixXcd::Identity(NS,NS)));  
  modes.resize(N_new, zeromode); 
  skip_list.resize(N_new, false);
}

void ModePropagatorGenerator_SingleModes::check_validity()const{
/*
  int N=get_N();
  if(HE.rows()!=N*rho_init.rows()){
    std::cerr<<"SingleMode: Dimension of HE must be factorizable into N_sys*N_env, whereas rho_E_init must be of dimension N_env!"<<std::endl;
    std::cerr<<"Hint: Please check parameters 'add_single_mode [HE] [rho_E_init]'!"<<std::endl;
    throw DummyException();
  }
*/
}

void ModePropagatorGenerator_SingleModes::add_single_mode(const std::vector<Eigen::MatrixXcd> & ops){
  if(ops.size()<2){
    std::cerr<<"ModePropagatorGenerator_SingleModes: ops.size()<2!"<<std::endl;
    throw DummyException();
  }
  Eigen::MatrixXcd HE=ops[0];
  Eigen::MatrixXcd rho_init=ops[1];
  std::vector<Eigen::MatrixXcd> env_ops;
  for(size_t i=2; i<ops.size(); i++)env_ops.push_back(ops[i]);
  modes.push_back(ModePropagatorPtr(new ModePropagator(rho_init, HE, env_ops )));
//  check_validity();
  if(skip_list.size()<1){skip_was_set=true;}
  skip_list.resize(modes.size(),false);
}

void ModePropagatorGenerator_SingleModes::add_single_mode(const std::vector<std::string> & str){
  std::vector<Eigen::MatrixXcd> ops;
  std::cout<<"Add single environment mode:"<<std::flush;
  for(size_t i=0; i<str.size(); i++){
    std::cout<<" '"<<str[i]<<"'"<<std::flush;
    ops.push_back(ReadExpression(str[i]));
  }
  std::cout<<std::endl;
  add_single_mode(ops);
}

void ModePropagatorGenerator_SingleModes::add_single_mode_from_file(const std::string &file, const std::vector<Eigen::MatrixXcd> & ops){
  if(ops.size()<1){
    std::cerr<<"ModePropagatorGenerator_SingleModes: ops.size()<1!"<<std::endl;
    throw DummyException();
  }
  Eigen::MatrixXcd rho_init=ops[0];
  std::vector<Eigen::MatrixXcd> env_ops;
  for(size_t i=1; i<ops.size(); i++)env_ops.push_back(ops[i]);

  Parameters param2;
  param2.add_from_file(file);
  modes.push_back(std::make_shared<ModePropagator>(FreePropagator(param2), rho_init, env_ops));
  if(skip_list.size()<1){skip_was_set=true;}
  skip_list.resize(modes.size(),false);
}


void ModePropagatorGenerator_SingleModes::add_single_mode(std::vector<std::shared_ptr<ModePropagator> > & list){
  for(size_t i=0; i<list.size(); i++){
   if(list[i]){modes.push_back(list[i]);}
  }
  if(skip_list.size()<1){skip_was_set=true;}
  skip_list.resize(modes.size(),false);
}

void ModePropagatorGenerator_SingleModes::add_single_mode_from_file(const std::vector<std::string> & str){
  if(str.size()<2){
    std::cerr<<"single_mode_from_file: needs at least 2 arguments: filename bath_init!"<<std::endl;
    throw DummyException();
  }
  std::vector<Eigen::MatrixXcd> ops;
  for(size_t i=1; i<str.size(); i++){
    std::cout<<" '"<<str[i]<<"'"<<std::flush;
    ops.push_back(ReadExpression(str[i]));
  }
  add_single_mode_from_file(str[0],ops);
}

EnvironmentOperators ModePropagatorGenerator_SingleModes::get_env_ops(int k) const{
  if(k>=modes.size()){
    std::cerr<<"ModePropagatorGenerator_SingleModes: get_env_ops(k) out of bounds!"<<std::endl;
    throw DummyException();
  }
  return EnvironmentOperators(modes[k]->env_ops);
}

Eigen::MatrixXcd ModePropagatorGenerator_SingleModes::get_bath_init(int k)const{
  if(k>=modes.size()){
    std::cerr<<"ModePropagatorGenerator_SingleModes: get_bath_init(k) out of bounds!"<<std::endl;
    throw DummyException();
  }
  return modes[k]->bath_init;
}
ModePropagatorPtr ModePropagatorGenerator_SingleModes::get_ModePropagator(int k)const{
    if(k>=modes.size()){
    std::cerr<<"ModePropagatorGenerator_SingleModes: get_bath_init(k) out of bounds!"<<std::endl;
    throw DummyException();
  }
  return modes[k];
}


}
