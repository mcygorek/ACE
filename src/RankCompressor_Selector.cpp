#include "RankCompressor_Selector.hpp"
#include "Smart_Ptr.h"
#include "RankCompressorList.hpp"
#include "Parameters.hpp"

namespace ACE{

RankCompressor_Ptr RankCompressor_Selector(Parameters &param, bool verbose){
 
  bool use_rrqr=param.get_as_bool("use_rrqr",false);
  double threshold=param.get_as_double("threshold", 0);
  double sum_threshold=param.get_as_double("sum_threshold", 0);
  int compress_maxk=param.get_as_size_t("compress_maxk", 0);
  double eps2=param.get_as_double("approxSVD",0);
  
  if(eps2>0.){
    if(verbose)std::cout<<"Using compression: approxSVD: "<<eps2<<" threshold: "<<threshold<<" sum_threshold: "<<sum_threshold<<" compress_maxk: "<<compress_maxk<<std::endl;
    return RankCompressor_Ptr(new RankCompressor_ApproxSVD(param));
  }else if(use_rrqr){
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

}//namespace
