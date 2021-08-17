#ifndef RANK_COMPRESSOR_SELECTOR_DEFINED_H
#define RANK_COMPRESSOR_SELECTOR_DEFINED_H

#include "Smart_Ptr.h"
#include "RankCompressor.h"
#include "Parameters.h"

typedef Smart_Ptr<RankCompressor> RankCompressor_Ptr;

RankCompressor_Ptr RankCompressor_Selector(Parameters &param, bool verbose=false){
 
  bool use_rrqr=param.get_as_bool("use_rrqr",false);
  double threshold=param.get_as_double("threshold", 0);
  double sum_threshold=param.get_as_double("sum_threshold", 0);
  int compress_maxk=param.get_as_size_t("compress_maxk", 0);
  
  if(use_rrqr){
    if(verbose)std::cout<<"Using compression: RRQR: threshold: "<<threshold<<std::endl;
    return RankCompressor_Ptr(new RankCompressor_RRQR(threshold));
#ifdef ALLOW_NEGATIVE_THRESHOLD
  }if(true){
#else
  }if(threshold>0 || compress_maxk>0 || sum_threshold>0){
#endif
    if(verbose)std::cout<<"Using compression: SVD: threshold: "<<threshold<<" sum_threshold: "<<sum_threshold<<" compress_maxk: "<<compress_maxk<<std::endl;
    return RankCompressor_Ptr(new RankCompressor_SVD(param));
  }else{
    if(verbose)std::cout<<"Using compression: None"<<std::endl;
    return RankCompressor_Ptr(new RankCompressor_None());
  }
}


#endif
