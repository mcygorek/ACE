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
  std::string outfile=param.get_as_string_check("outfile");
  int n_extract=param.get_as_int("n", -1);
  double threshold=param.get_as_double("threshold",0);
  bool print_timesteps=param.get_as_bool("print_timesteps", false);
  bool skip_sweep_backward=param.get_as_bool("skip_sweep_backward",false);

  int final_sweep_n=param.get_as_int("final_sweep_n",0);
  RankCompressor_SVD compr_final; 
  {Parameters param2; param2.add_from_prefix("final_sweep",param);
   compr_final.setup(param2);}

  InfluenceFunctional_OD IF(read_PT);
  
  if(IF.a.size()<1){
    std::cerr<<"IF.a.size()<1!"<<std::endl;
    exit(1);
  }
  if(n_extract<0)n_extract=IF.a.size()/2;
  if(n_extract>=IF.a.size()-1){
    std::cerr<<"n_extract>=IF.a.size()-1 ("<<n_extract<<" vs. "<<IF.a.size()<<")!"<<std::endl;
    exit(1);
  }

IF.check_env_dims();
  for(int loop=0; loop<final_sweep_n; loop++){
//    std::cout<<"Loop: "; compr_final.print_info(); std::cout<<std::endl;
    std::cout<<"Loop "<<loop<<"/"<<final_sweep_n<<": Sweeping backward"<<std::endl;
    sweep_high_to_low(IF, compr_final, print_timesteps);

    std::cout<<"Loop "<<loop<<"/"<<final_sweep_n<<": Sweeping forward"<<std::endl;
    sweep_low_to_high(IF, compr_final, print_timesteps);
  }
IF.check_env_dims();

  if(!skip_sweep_backward){
    std::cout<<"Sweeping backward"<<std::endl;
    RankCompressor_SVD compr(param);
    sweep_high_to_low(IF, compr, print_timesteps);
  }

  std::cout<<"To forward normal form and ";
  std::cout<<"extract n="<<n_extract<<"/"<<IF.a.size()<<std::endl;
  std::ofstream ofs(outfile.c_str());

  std::vector<Eigen::VectorXd> forwardNF=make_forwardNF(IF, threshold, print_timesteps);

  for(int r=0; r<forwardNF[n_extract].rows(); r++){
    ofs<<forwardNF[n_extract](r)<<std::endl;
  }

  return 0;
}


