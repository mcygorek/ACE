#include "ModePropagatorGenerator.hpp"
#include "ModePropagatorGenerator_Interweave.hpp"
#include "MPG_Selector.hpp"
#include "Parameters.hpp"

namespace ACE{
  std::string ModePropagatorGenerator_Interweave::name()const{
    return std::string("Interweave");
  } 

  void ModePropagatorGenerator_Interweave::setup(const std::vector<std::string> & str, Parameters &param){

    if(str.size()<1){
      std::cerr<<"ModePropagatorGenerator_Interweave: setup needs at least 1 filename"<<std::endl;
      exit(1);
    }

    //make mpgs list:
    mpgs.clear();
    for(size_t i=0; i<str.size(); i++){
      Parameters param(str[i]);
      std::vector<std::shared_ptr<ModePropagatorGenerator> >  this_mpgs = MPG_Selector(param);
      mpgs.insert(mpgs.end(), this_mpgs.begin(), this_mpgs.end());
    }

    //generate "which" list:
    which.clear();
    std::string order=param.get_as_string("interweave_order");
    if(order=="sequential"){
      for(int i=0; i<(int)mpgs.size(); i++){
        for(int j=0; j<mpgs[i]->get_N_modes(); j++){
          which.push_back({i, j});
        }
      }
    }else{
      int max_modes=0; 
      for(size_t i=0; i<mpgs.size(); i++){
        if(mpgs[i]->get_N_modes() > max_modes)max_modes=mpgs[i]->get_N_modes();
      }
      for(int j=0; j<max_modes; j++){
        for(int i=0; i<(int)mpgs.size(); i++){
          if(j<mpgs[i]->get_N_modes()){
            which.push_back({i, j});
          }
        }
      }
    }

    std::cout<<"ModePropagatorGenerator_Interweave: setup:"<<std::endl;
    print_which();
    setup_skip(param);
  }
 
  void ModePropagatorGenerator_Interweave::check_exists(int k, const std::string &context)const{

    if(k<0||k>=get_N_modes()){
      std::cerr<<"ModePropagatorGenerator_Interweave: k<0||k>=get_N_modes() ("<<k<<" vs. "<<get_N_modes()<<")!"<<std::endl; 
      throw DummyException();
    }

    if(which[k].first<0 || which[k].first>=mpgs.size()){
      std::cerr<<"ModePropagatorGenerator_Interweave: which[k].first>=mpgs.size() ("<<which[k].first<<" vs. "<<mpgs.size()<<")!"<<std::endl; 
      throw DummyException();
    }
  }

  void ModePropagatorGenerator_Interweave::print_which(std::ostream &os)const{
    os<<"which.size()="<<which.size()<<std::endl;
    for(int k=0; k<which.size(); k++){
      os<<k<<": ("<<which[k].first<<","<<which[k].second<<") ";
      if(mpgs.size()>=which[k].first){ os<<mpgs[which[k].first]->name(); }
      os<<std::endl;
    }
  }


  int ModePropagatorGenerator_Interweave::get_N()const{
    if(mpgs.size()<1)return 0;
    return mpgs[0]->get_N(); 
  }

  std::vector<Eigen::MatrixXcd> ModePropagatorGenerator_Interweave::get_env_ops(int k)const{
    check_exists(k,"get_env_ops");
    return mpgs[which[k].first]->get_env_ops(which[k].second);
  }

  Eigen::MatrixXcd ModePropagatorGenerator_Interweave::get_bath_init(int k)const{
    check_exists(k,"get_bath_init");
    return mpgs[which[k].first]->get_bath_init(which[k].second);
  }

  int ModePropagatorGenerator_Interweave::get_mode_dim(int k)const{
    check_exists(k,"get_mode_dim");
    return mpgs[which[k].first]->get_mode_dim(which[k].second);
  }

  ModePropagatorPtr ModePropagatorGenerator_Interweave::getModePropagator(int k)const{
    check_exists(k,"getModePropagator");
    return mpgs[which[k].first]->getModePropagator(which[k].second);
  }

}//namespace
