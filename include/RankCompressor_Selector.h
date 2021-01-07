#ifndef RANK_COMPRESSOR_SELECTOR_DEFINED_H
#define RANK_COMPRESSOR_SELECTOR_DEFINED_H

#include "Smart_Ptr.h"
#include "RankCompressor.h"
#include "Parameters.h"

typedef Smart_Ptr<RankCompressor> RankCompressor_Ptr;

RankCompressor_Ptr RankCompressor_Selector(Parameters &param){
 
  bool use_rrqr=param.get_as_bool("use_rrqr",false);
  double threshold=param.get_as_double("threshold", 0);
  int compress_maxk=param.get_as_size_t("compress_maxk", 0);
  bool reortho=param.get_as_bool("reorthogonalize",false);
  
  if(use_rrqr){
    std::cout<<"Using compression: RRQR: threshold: "<<threshold<<std::endl;
    return RankCompressor_Ptr(new RankCompressor_RRQR(threshold));
#ifdef ALLOW_NEGATIVE_THRESHOLD
  }if(true){
#else
  }if(threshold>0 || compress_maxk>0){
#endif
    std::cout<<"Using compression: SVD: threshold: "<<threshold<<" compress_maxk: "<<compress_maxk<<std::endl;
    return RankCompressor_Ptr(new RankCompressor_SVD(threshold,compress_maxk,reortho));
  }else{
    std::cout<<"Using compression: None"<<std::endl;
    return RankCompressor_Ptr(new RankCompressor_None());
  }
}


#endif
