#ifndef IF_OD_PLAN_DEFINED_H
#define IF_OD_PLAN_DEFINED_H

#include "InfluenceFunctional_OD.hpp"
#include "Closure_Ops.hpp"
#include "Env_State_Filter.hpp"
//#include "ProcessTensor_real.hpp"
#include "RankCompressor_Selector.hpp"
#include "MPG_Selector.hpp"
#include "RankCompressor_Selector.hpp"
//#include "ModePropagatorGeneratorList.hpp"

/* 
   Read from Parameter what process tensors are to be calculated.
   Later, the Plan can be executed.
*/
namespace ACE{
class Parameters;

class IF_OD_Plan{
public:
  TimeGrid tgrid;
  bool use_dict;
  double dict_zero;
  int n_coarse;
  int factorization;
  bool compress_trafo_use_ortho;

  bool IF_print_timesteps;
  int intermediate_sweeps;
  RankCompressor_SVD_real intermediate_sweeps_compressor;
  int final_sweeps; //NOTE: that's for the ProcessTensor_real, not the complex IF_OD! (cf. "final_sweep" without trailing "s")
  RankCompressor_SVD_real final_sweeps_compressor;

  RankCompressor_Ptr compressor;

  bool use_realPT;
  RankCompressor_SVD_real compressor_real;
  Closure_Ops closure_ops;

  int final_sweep_n; //NOTE: that's for the complex IF_OD, not the ProcessTensor_real! (cf. "final_sweeps" with a trailing "s")
  RankCompressor_Ptr final_sweep_compressor;

  bool use_Gaussian;
  int n_extra;
  DiagBB Gaussian_DiagBB;

  int dim;
  std::string write_PT;
  ReadPT_struct read_PT;
  std::vector<ReadPT_struct> multi_PT;

  std::vector<std::shared_ptr<ModePropagatorGenerator> > mpgs;
  std::vector<std::shared_ptr<FreePropagator> > PT_apply_prop;

  std::shared_ptr<Env_State_Filter> env_state_filter;

  std::string print_dims_to_file;

  std::vector<std::pair<int, std::string> > trafo_chain_files;

  void check_consistency(int check_sys_dim=0)const;

  void setup(Parameters &param, int check_sys_dim=0);

  //===================================================

  inline bool is_trivial()const{
    if(!!read_PT)return false;
    if(use_Gaussian)return false;
    if(mpgs.size()>0)return false;
    return true;
  }

  std::vector<std::shared_ptr<InfluenceFunctional_OD> > execute();

  IF_OD_Plan(Parameters &param, int check_sys_dim=0);
  IF_OD_Plan();
};


}//namespace
#endif
