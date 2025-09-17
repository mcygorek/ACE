#include "TruncatedSVD.hpp"
#include "TruncationLayout.hpp"
#include <Eigen/SVD>
#include "DummyException.hpp"
#include "ReadTable.hpp"

namespace ACE{

template <typename T> void TruncationLayout_T<T>::print_info(std::ostream &ofs)const{
  ofs<<"base_threshold="<<base_threshold;
  ofs<<" base_maxk="<<base_maxk;
  ofs<<" base_mink="<<base_mink;
  ofs<<" keep="<<keep;
  ofs<<" forward_threshold_ratios="<<forward_threshold_ratio;
  ofs<<" backward_threshold_ratios="<<backward_threshold_ratio;
  ofs<<" select_threshold_ratios="<<select_threshold_ratio;
  ofs<<" threshold_range_factor="<<threshold_range_factor;
  ofs<<" QR_after_combine=";if(use_QR)ofs<<"true";else ofs<<"false";
  if(base_Tikhonov>0.)ofs<<" Tikhonov="<<base_Tikhonov;
  ofs<<" intermediate_sweep_n="<<intermediate_sweep_n;
  ofs<<" final_sweep_n=";
  if(final_sweep_half){ ofs<<final_sweep_n+0.5; }else{ ofs<<final_sweep_n; }
  ofs<<" final_sweep_threshold="<<final_sweep_threshold;
  ofs<<" final_sweep_maxk="<<final_sweep_maxk;
}


template <typename T> void TruncationLayout_T<T>::setup(Parameters &param){
  base_threshold = param.get_as_double("threshold",0.);
  base_maxk = param.get_as_int("compress_maxk",0);
  base_mink = param.get_as_int("compress_mink",0);
  keep = param.get_as_double("compress_keep",-1);
  base_Tikhonov = param.get_as_double("compress_Tikhonov",0.);

  use_QR = param.get_as_bool("QR_after_combine_twice",false);
  forward_threshold_ratio = param.get_as_double("forward_threshold_ratio",1.);
  backward_threshold_ratio = param.get_as_double("backward_threshold_ratio",1.);
  select_threshold_ratio = param.get_as_double("select_threshold_ratio",1.);
  threshold_range_factor = param.get_as_double("threshold_range_factor",1.);


  intermediate_sweep_n = param.get_as_size_t("intermediate_sweep_n",0);
  double fin_n= param.get_as_double("final_sweep_n", intermediate_sweep_n);
  final_sweep_n = fin_n;
  final_sweep_half = (fin_n-final_sweep_n>0.49);
  //manually override default when use_Gaussian_repeat is set:
  if(param.get_as_bool("use_Gaussian_repeat") && !param.is_specified("final_sweep_n")){final_sweep_half=true;}

  final_sweep_threshold = param.get_as_double("final_sweep_threshold",base_threshold);
  final_sweep_maxk = param.get_as_int("final_sweep_compress_maxk",base_maxk);
}

template class TruncationLayout_T<std::complex<double> >;
template class TruncationLayout_T<double>;

}//namespace
