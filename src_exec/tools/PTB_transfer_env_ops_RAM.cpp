#include "InfluenceFunctional_OD.hpp"
#include "Parameters.hpp"
#include "BlockCombine_Basics.hpp"
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <vector>

using namespace ACE;

int main(int args, char** argv){
  Parameters param(args, argv);

  std::string read_PT=param.get_as_string_check("read_PT");
  std::string multi_PT=param.get_as_string_check("multi_PT");
  std::string outfile=param.get_as_string("outfile");
  std::string write_PT=param.get_as_string("write_PT");
  int n_extract=param.get_as_int("n", -1);
  double threshold=param.get_as_double("threshold",0);
  bool print_timesteps=param.get_as_bool("print_timesteps", false);
  bool skip_sweep_backward=param.get_as_bool("skip_sweep_backward",false);
  int fix_phase=param.get_as_size_t("fix_phase",1);

  int final_sweep_n=param.get_as_int("final_sweep_n",0);
  double final_sweep_threshold=param.get_as_double("final_sweep_threshold",threshold);

  InfluenceFunctional_OD IF(read_PT);
  InfluenceFunctional_OD IF2(multi_PT);
  IF.check_env_dims();
  IF2.check_env_dims();

  
  if(IF.a.size()!=IF2.a.size()){
    std::cerr<<"IF.a.size()!=IF2.a.size() ("<<IF.a.size()<<" vs. "<<IF2.a.size()<<")!"<<std::endl;
    exit(1);
  }
  if(IF.a.size()<1){
    std::cerr<<"IF.a.size()<1!"<<std::endl;
    exit(1);
  }
  if(n_extract<0)n_extract=IF.a.size()/2;
  if(n_extract>=IF.a.size()-1){
    std::cerr<<"n_extract>=IF.a.size()-1 ("<<n_extract<<" vs. "<<IF.a.size()<<")!"<<std::endl;
    exit(1);
  }

  if(!skip_sweep_backward){
    std::cout<<"Sweeping backward"<<std::endl;
    RankCompressor_SVD compr(threshold);
//    IF.check_env_dims();
//    IF2.check_env_dims();
    sweep_high_to_low(IF, compr, print_timesteps);
    sweep_high_to_low(IF2, compr, print_timesteps);
  }

  for(int loop=0; loop<final_sweep_n; loop++){
    RankCompressor_SVD compr(threshold);
    std::cout<<"Loop "<<loop<<"/"<<final_sweep_n<<": Sweeping forward"<<std::endl;
    sweep_low_to_high(IF, compr, print_timesteps);
    sweep_low_to_high(IF2, compr, print_timesteps);

    std::cout<<"Loop "<<loop<<"/"<<final_sweep_n<<": Sweeping backward"<<std::endl;
    sweep_high_to_low(IF, compr, print_timesteps);
    sweep_high_to_low(IF2, compr, print_timesteps);
  }

  IF.calculate_closures();
  IF2.calculate_closures();

  std::cout<<"To forward normal form and ";
  std::cout<<"extract n="<<n_extract<<"/"<<IF.a.size()<<std::endl;
  std::cout<<"using fix_phase="<<fix_phase<<std::endl;

  std::vector<Eigen::VectorXd> forwardNF=make_forwardNF(IF, final_sweep_threshold, print_timesteps, fix_phase);

  std::vector<Eigen::VectorXd> forwardNF2=make_forwardNF(IF2, final_sweep_threshold, print_timesteps, fix_phase);


  for(int r=0; r<forwardNF[n_extract].rows(); r++){
    std::cout<<forwardNF[n_extract](r)<<" ";
  }
  std::cout<<std::endl;


  for(int r=0; r<forwardNF2[n_extract].rows(); r++){
    std::cout<<forwardNF2[n_extract](r)<<" ";
  }
  std::cout<<std::endl;

  if(outfile!=""){
    std::ofstream ofs(outfile.c_str());
    for(int r=0; r<forwardNF[n_extract].rows() || r<forwardNF2[n_extract].rows(); r++){
      if(r<forwardNF[n_extract].rows()){
        ofs<<forwardNF[n_extract](r)<<" ";
      }else{
        ofs<<"nan ";
      }
      if(r<forwardNF2[n_extract].rows()){
        ofs<<forwardNF2[n_extract](r)<<" ";
      }else{
        ofs<<"nan ";
      }
      ofs<<std::endl;
    }
  }

  if(write_PT!=""){
    if(IF.env_ops.size()<1){
      std::cerr<<"IF.env_ops.size()<1!"<<std::endl;
      exit(1);
    }
    IF2.env_ops=std::vector<std::vector<Eigen::VectorXcd> >(IF2.a.size());
    for(int n=0; n<IF.a.size(); n++){
      IF2.env_ops[n]=std::vector<Eigen::VectorXcd>(IF.env_ops[n].size());
      for(int o=0; o<IF.env_ops[n].size(); o++){
        IF2.env_ops[n][o]=Eigen::VectorXcd::Zero(IF2.a[n].dim_d2);
        for(int r=0; r<IF.env_ops[n][o].rows() && r<IF2.env_ops[n][o].rows(); r++){
          IF2.env_ops[n][o](r)=IF.env_ops[n][o](r);
        }
      }
    }
    IF2.write_binary(write_PT);
  }

  return 0;
}


