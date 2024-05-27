#include "InfluenceFunctional_Repeat.hpp"
#include "IF_OD_Abstract.hpp"
#include "ModePropagatorGenerator.hpp"
#include "ProcessTensor_real.hpp"
#include "BinaryReader.hpp"
#include "MPS_Matrix.hpp"

namespace ACE{

  void InfluenceFunctional_Repeat::print_debug(std::ostream &ofs)const{
    ofs<<"InfluenceFunctional_Repeat: a.size()="<<a.size();
    ofs<<" c.size()="<<c.size()<<" env_ops.size()="<<env_ops.size()<<":";
    for(size_t i=0; i<env_ops.size(); i++)std::cout<<" "<<env_ops[i].size();
    std::cout<<std::endl;
  }

  //Implementation of IF_OD_Abstract
  const MPS_Matrix & InfluenceFunctional_Repeat::get_a(int n)const{
    if(n<0){
      std::cerr<<"InfluenceFunctional_OD: get_a out of bounds "<<n<<"/"<<a.size()<<std::endl;
      exit(1);
    }
    if(n==0){
      return a[0];
    }else{
      return a[1];
    }
  }

  const Eigen::VectorXcd & InfluenceFunctional_Repeat::get_c(int n)const{
    if(n<0){
      std::cerr<<"InfluenceFunctional_OD: get_c out of bounds "<<n<<"/"<<c.size()<<std::endl;
      exit(1);
    }
    if(n==0){
      return c[0];
    }else{
      return c[1];
    }
  }
  const std::vector<Eigen::VectorXcd> & InfluenceFunctional_Repeat::get_env_ops(int n)const{
    if(n<0){
      std::cerr<<"InfluenceFunctional_OD: Trying to access env_ops out of bounds "<<n<<"/"<<env_ops.size()<<std::endl;
      exit(1);
    }
    if(n==0){ 
      return env_ops[0];
    }else{
      return env_ops[1];
    }
  }

  void InfluenceFunctional_Repeat::check_within_limits(int n)const{
    if(n<0){
        std::cerr<<"InfluenceFunctional_Repeat::check_within_limits: n<0!"<<std::endl;
        exit(1);
    }
  }
  

  //Genuine IF_OD functions:
  void InfluenceFunctional_Repeat::calculate_closures(){
    if(a.size()!=3){
      std::cerr<<"InfluenceFunctional_Repeat::calculate_closures: a!=3!"<<std::endl;
      exit(1);
    }

    int N=dict.get_N();
//    int NL=dict.get_NL();

    c.resize(a.size());
    c.back().resize(1);
    c.back()(0)=1;

    for(int n=(int)c.size()-2; n>=0; n--){
      c[n]=Eigen::VectorXcd::Zero(a[n].dim_d2);
      for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
          int i_ind=dict(((i*N+i)*N+j)*N+j);
          if(i_ind<0)continue;
          for(int d1=0; d1<a[n+1].dim_d1; d1++){
            for(int d2=0; d2<a[n+1].dim_d2; d2++){
              c[n](d1)+=c[n+1](d2)*a[n+1](i_ind, d1, d2)/((double)N);
            }
          }
        }
      }
    }
  }

  void InfluenceFunctional_Repeat::read_binary(const std::string &fname){
    std::ifstream ifs(fname.c_str());
    if(!ifs.good()){
      std::cerr<<"Cannot read InfluenceFunctional_Repeat file '"<<fname<<"'!"<<std::endl;
      exit(1);
    }
    if(binary_read_fixedSizeString(ifs, 4, "InfluenceFunctional_Repeat")!="RePT"){
      std::cerr<<"Cannot interpret '"<<fname<<"' as InfluenceFunctional_Repeat file!"<<std::endl;
      exit(1);
    }
    dict.read_binary(ifs);
    
    c.resize(binary_read_int(ifs, "IF_Repeat"));
    for(size_t i=0; i<c.size(); i++){
      c[i]=binary_read_EigenMatrixXcd(ifs, "IF_Repeat");
    }
  
    env_ops.resize(binary_read_int(ifs, "IF_Repeat"));
    for(size_t i=0; i<env_ops.size(); i++){
      env_ops[i].resize(binary_read_int(ifs, "IF_Repeat"));
      for(size_t j=0; j<env_ops[i].size(); j++){
        env_ops[i][j]=binary_read_EigenMatrixXcd(ifs, "IF_Repeat");
      }
    }
      
    MPS::read_binary(ifs); 
    
    calculate_closures();
    print_debug();
  }

  void InfluenceFunctional_Repeat::write_binary(const std::string &fname)const{
    std::ofstream ofs(fname.c_str());
    binary_write_fixedSizeString(ofs,4,"RePT");
    dict.write_binary(ofs);
   
    binary_write_int(ofs,c.size());
    for(size_t i=0; i<c.size(); i++){
      binary_write_EigenMatrixXcd(ofs, c[i]);
    }
  
    binary_write_int(ofs,env_ops.size());
    for(size_t i=0; i<env_ops.size(); i++){
      binary_write_int(ofs,env_ops[i].size());
      for(size_t j=0; j<env_ops[i].size(); j++){
        binary_write_EigenMatrixXcd(ofs, env_ops[i][j]);
      }
    }
   
    MPS::write_binary(ofs);
  }

  InfluenceFunctional_Repeat::InfluenceFunctional_Repeat(int N){
    dict.set_default(N);
    a.resize(3);
    c.resize(3);
    env_ops.resize(3);
  }
  InfluenceFunctional_Repeat::InfluenceFunctional_Repeat(const std::string &fname){
    read_binary(fname);
  }

}//namespace
