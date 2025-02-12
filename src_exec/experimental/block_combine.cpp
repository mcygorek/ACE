#include "PCH.hpp"
#include "ACE.hpp"
#include "BlockCombine_Basics.hpp"
#include "BlockCombine.hpp"

using namespace ACE;

int main(int args, char** argv){

  Parameters param(args, argv, true);

  TimeGrid tgrid(param);
  std::string read_PT = param.get_as_string("read_PT");
  std::vector<std::string> multi_PT = param.get_all_strings("multi_PT");
  
  std::string write_PT = param.get_as_string_check("write_PT");
  double threshold = param.get_as_double_check("threshold");
  double dict_zero = param.get_as_double("dict_zero",0);

  bool verbose=param.get_as_bool("verbose",false);
  bool IF_print_timesteps=param.get_as_bool("IF_print_timesteps",false);
  RankCompressor_Ptr compressor=RankCompressor_Selector(param, true);

//  int intermediate_sweep_n=param.get_as_int("intermediate_sweep_n", 0);

  bool noSVD=param.get_as_bool("noSVD",false);
  bool ignore_multi_env_ops=param.get_as_bool("ignore_multi_env_ops",false);

  int final_sweep_n=param.get_as_size_t("final_sweep_n",0);
  double final_sweep_threshold=param.get_as_double("final_sweep_threshold",threshold);

  int default_dim=2;
  std::vector<std::shared_ptr<ModePropagatorGenerator> > mpgs=MPG_Selector(param);
  if(mpgs.size()>0){
    default_dim=mpgs[0]->get_N();
  }

  InfluenceFunctional_OD IF(tgrid, default_dim);
  if(read_PT!=""){
    IF.read_binary(read_PT);
  }
  
  BlockCombine BC(param);
  BC.verbose=true;

  for(auto & mpg : mpgs){
    BC.BlockCombine_MPG(IF, *mpg.get());
  }  

  for(const std::string &str : multi_PT){
    InfluenceFunctional_OD IF2(str);
    if(ignore_multi_env_ops){
      IF2.env_ops.clear();
    }
//    Block_Combine_and_Sweep_Back(IF, IF2, threshold, verbose, IF_print_timesteps);
    Block_Combine_and_Sweep_NF(IF, IF2, threshold, verbose, IF_print_timesteps, noSVD);
  }

  for(int i=0; i<final_sweep_n; i++){
    double thr=threshold;
    if(i==final_sweep_n-1)thr=final_sweep_threshold;
    RankCompressor_SVD compr(thr);

    std::cout<<"Final sweep "<<i<<"/"<<final_sweep_n<<" forward (thr="<<thr<<")"<<std::endl;
    sweep_low_to_high(IF, compr, IF_print_timesteps);
    if(verbose){
      std::cout<<"After compression: IF.get_max_dim()="<<IF.get_max_dim()<<std::endl;
    }

    std::cout<<"Final sweep "<<i<<"/"<<final_sweep_n<<" backward (thr="<<thr<<")"<<std::endl;
    sweep_high_to_low(IF, compr, IF_print_timesteps);
    if(verbose){
      std::cout<<"After compression: IF.get_max_dim()="<<IF.get_max_dim()<<std::endl;
    }
  }  

  IF.calculate_closures();
  IF.write_binary(write_PT);

  return 0;
}

