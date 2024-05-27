#include "BlockCombine_Basics.hpp"
#include "BlockCombine.hpp"
#include "MPG_Selector.hpp"

namespace ACE{

void BlockCombine::BlockCombine_MPG_sequential(InfluenceFunctional_OD &IF, 
                                    ModePropagatorGenerator &mpg) const{

  RankCompressor_SVD compressor(threshold);

  for(int k=mpg.first(); k<mpg.get_N_modes(); k=mpg.next(k)){
    if(verbose){std::cout<<"Mode: k="<<k<<std::endl;}
    ModePropagatorPtr mpp=mpg.getModePropagator(k);

    InfluenceFunctional_OD IF2(*mpp.get(), tgrid, dict_zero);
 
//    Block_Combine_and_Sweep_Back(IF, IF2, threshold, verbose, IF_print_timesteps);
    Block_Combine_and_Sweep_NF(IF, IF2, threshold, verbose, IF_print_timesteps);
  }
}


// Divide and conquer strategy:
// decompose into hierarchy:
// lvl 0:   0 1 2 3 4 5 6 7 
// lvl 1:    x   x   x   x  
// lvl 2:      x       x
// lvl 3:          x
// solve by recursion:

InfluenceFunctional_OD BlockCombine::BlockCombine_MPG_tree_get(
       int level, int first_elem, ModePropagatorGenerator &mpg) const{


  if(level==0){
    if(first_elem>=mpg.get_N_modes()){
      std::cerr<<"BlockCombine_MPG_tree_get: level="<<level<<" first_elem="<<first_elem<<": first_element>=mpg.get_N_modes()="<<mpg.get_N_modes()<<"!"<<std::endl;
      exit(1);
    }

    if(mpg.skip_list[first_elem]){
      return InfluenceFunctional_OD(tgrid, mpg.get_N());
    }else{
      if(verbose){
        std::cout<<"level: "<<level<<" first_elem: "<<first_elem<<std::endl;
      }
      ModePropagatorPtr mpp=mpg.getModePropagator(first_elem);
      return InfluenceFunctional_OD(*mpp.get(), tgrid, dict_zero);
    }

  }else{
    InfluenceFunctional_OD IF = 
        BlockCombine_MPG_tree_get(level-1, 2*first_elem, mpg);

    InfluenceFunctional_OD IF2 = 
        BlockCombine_MPG_tree_get(level-1, 2*first_elem+1, mpg);

    if(verbose){
      std::cout<<"level: "<<level<<" first_elem: "<<first_elem<<std::endl;
    }
//    Block_Combine_and_Sweep_Back(IF, IF2, threshold, verbose, IF_print_timesteps);

    Block_Combine_and_Sweep_NF(IF, IF2, threshold, verbose, IF_print_timesteps);

    RankCompressor_SVD compressor(threshold);
    for(int im=0; im<intermediate; im++){
      sweep_low_to_high(IF, compressor, IF_print_timesteps);
      if(verbose)std::cout<<"intermediate (forward): "<<im<<"/"<<intermediate<<std::flush;
      if(verbose)std::cout<<" max_dim="<<IF.get_max_dim()<<std::endl;

      sweep_high_to_low(IF, compressor, IF_print_timesteps);
      if(verbose)std::cout<<"intermediate (backward): "<<im<<"/"<<intermediate<<std::flush;
      if(verbose)std::cout<<" max_dim="<<IF.get_max_dim()<<std::endl;
    }   

    return IF;
  }
}

void BlockCombine::BlockCombine_MPG_tree(InfluenceFunctional_OD &IF, 
                                         ModePropagatorGenerator &mpg)const{

  RankCompressor_SVD compressor(threshold);
  int N_modes=mpg.get_N_modes();
  if(N_modes<1)return;
 
  int N_hierarchy=1;
  {int N_shift=N_modes;
    while(N_shift>1){
      N_hierarchy++;
      N_shift=N_shift>>1;
    }
  }
  if(N_modes!=pow(2, N_hierarchy-1)){
    std::cerr<<"BlockCombine_MPG_tree: N_modes="<<N_modes<<" != 2^"<<N_hierarchy-1<<"="<<pow(2, N_hierarchy-1)<<" => not a power of 2!"<<std::endl;
    exit(1);
  }
  std::cout<<"N_modes="<<N_modes<<"=2^"<<N_hierarchy-1<<std::endl;

  IF=BlockCombine_MPG_tree_get(N_hierarchy-1, 0, mpg);

}

void BlockCombine::setup(Parameters &param){
  tgrid = TimeGrid(param);
  threshold = param.get_as_double("threshold");
  intermediate = param.get_as_int("sweep_intermediate_n", 0);
  dict_zero = param.get_as_double("dict_zero",0);
  verbose = false;
  IF_print_timesteps = param.get_as_bool("IF_print_timesteps",false);
  use_BC_tree = param.get_as_bool("use_BC_tree", false);
}

}//namespace
