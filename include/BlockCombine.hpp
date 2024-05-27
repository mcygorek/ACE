#ifndef ACE_BLOCKCOMBINE_MPG_DEFINED_H
#define ACE_BLOCKCOMBINE_MPG_DEFINED_H

#include "InfluenceFunctional_OD.hpp"
#include "ModePropagatorGenerator.hpp"
#include "Parameters.hpp"

namespace ACE{

class BlockCombine{
public:
  TimeGrid tgrid;

  double threshold;
  int intermediate;
  double dict_zero; 
  bool verbose;
  bool IF_print_timesteps;
  bool use_BC_tree;


  void BlockCombine_MPG_sequential(InfluenceFunctional_OD &IF, 
                                   ModePropagatorGenerator &mpg) const;


  InfluenceFunctional_OD BlockCombine_MPG_tree_get(
       int level, int first_elem, ModePropagatorGenerator &mpg) const;

  void BlockCombine_MPG_tree(InfluenceFunctional_OD &IF, 
                             ModePropagatorGenerator &mpg) const;

  inline void BlockCombine_MPG(InfluenceFunctional_OD &IF, 
                               ModePropagatorGenerator &mpg) const{
    if(use_BC_tree){
      BlockCombine_MPG_tree(IF, mpg);
    }else{
      BlockCombine_MPG_sequential(IF, mpg);
    }
  }

  void setup(Parameters &param);

  inline BlockCombine(Parameters &param){ 
    setup(param);
  }
  inline BlockCombine(){
    Parameters param;
    setup(param);
  }
  
};
}//namespace
#endif
