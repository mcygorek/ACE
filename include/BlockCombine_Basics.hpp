#ifndef ACE_BLOCKCOMBINE_BASICS_DEFINED_H
#define ACE_BLOCKCOMBINE_BASICS_DEFINED_H

#include "InfluenceFunctional_OD.hpp"
#include "RankCompressor.hpp"

namespace ACE{

void BlockCombine_FixedOrder(InfluenceFunctional_OD & IF1, const InfluenceFunctional_OD & IF2, double threshold, bool print_timesteps=false);

void BlockCombine_AlternateOrder(InfluenceFunctional_OD & IF1, const InfluenceFunctional_OD & IF2, double threshold, bool print_timesteps=false);


//assumes forward normal form of IF1 and IF2; 
//compresses in forward direction after combination 
//=> Turns out not to be useful
void BlockCombineCompress_NF(
  InfluenceFunctional_OD & IF1, const std::vector<Eigen::VectorXd> &NF1, 
  const InfluenceFunctional_OD & IF2, const std::vector<Eigen::VectorXd> &NF2,  
  double threshold, bool print_timesteps=false, int *maxdim_intermediate=NULL);


//assumes forward normal form of IF1 and IF2; 
//compresses in backward direction after combination 
void BlockCombineCompressBack_NF(
  InfluenceFunctional_OD & IF1, const std::vector<Eigen::VectorXd> &NF1, 
  const InfluenceFunctional_OD & IF2, const std::vector<Eigen::VectorXd> &NF2,  
  double threshold, bool print_timesteps=false, int *maxdim_intermediate=NULL,
  bool noSVD=false);



void sweep_high_to_low(InfluenceFunctional_OD &IF, RankCompressor &compressor, bool print_timesteps=false);

void sweep_low_to_high(InfluenceFunctional_OD &IF, RankCompressor &compressor, bool print_timesteps=false);


//fix_phase: 0: don't fix it
//           1: phase of element with max. abs
//           2: phase of first element with max. abs >1e-6
std::vector<Eigen::VectorXd> make_forwardNF(InfluenceFunctional_OD &IF, double threshold, bool print_timesteps=false, int fix_phase=0);

void Block_Combine_and_Sweep_Back(InfluenceFunctional_OD &IF,
                                  InfluenceFunctional_OD &IF2,
                                  double threshold, bool verbose,
                                  bool IF_print_timesteps=false);

//brings IF1 and IF2 to forward normal form and combines them 
void Block_Combine_and_Sweep_NF(InfluenceFunctional_OD &IF,
                                InfluenceFunctional_OD &IF2,
                                double threshold, bool verbose,
                                bool IF_print_timesteps=false, 
                                bool noSVD=false);
}//namespace


#endif
