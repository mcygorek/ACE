#include "EDM_Filter.hpp"

namespace ACE{

bool EDM_Filter::passes( const std::complex<double> &val) const{
  if(std::abs(val)>threshold){
    return true;  
  }
  return false;
}
bool EDM_Filter::candidate_passes( const std::complex<double> &val) const{
  if(std::abs(val)>candidate_threshold){
    return true;  
  }
  return passes(val);
}
bool EDM_Filter::increment_passes( const std::complex<double> &val){
  double v=std::abs(val);
  if(v>max_in_round){
    max_in_round=v;
  }
  if(threshold_in_round>0. && v<threshold_in_round*max_in_round){
    return false;  
  }
  if(increment_threshold>0. && v<increment_threshold){
    return false;  
  }
  return true;
}
bool EDM_Filter::sparse_prop_passes( const std::complex<double> &val) const{
  if(std::abs(val)>sparse_prop_threshold){
    return true;  
  }
  return false;
}

void EDM_Filter::print_info(std::ostream &ofs)const{
  ofs<<"threshold "<<threshold;
  ofs<<", candidate_threshold "<<candidate_threshold;
  ofs<<", increment_threshold "<<increment_threshold;
  ofs<<", threshold_in_round "<<threshold_in_round;
  ofs<<", sparse_prop_threshold "<<sparse_prop_threshold;
  ofs<<", max_coeffs "<<max_coeffs;
  ofs<<std::endl;
}

void EDM_Filter::setup(Parameters &param){
  threshold=param.get_as_double("threshold");
  candidate_threshold=param.get_as_double("candidate_threshold", threshold);
  sparse_prop_threshold=param.get_as_double("sparse_prop_threshold");
  increment_threshold=param.get_as_double("increment_threshold");
  threshold_in_round=param.get_as_double("threshold_in_round");

  max_coeffs=param.get_as_int("max_coeffs");
  mc_max=param.get_as_double("mc_max", 0.1); 
  mc_min=param.get_as_double("mc_min", 1e-20); 
  mc_Nsteps=param.get_as_size_t("mc_Nsteps",200);
}
}//namespace
