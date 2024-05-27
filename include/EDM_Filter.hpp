#ifndef ACE_EDM_FILTER_DEFINED_H
#define ACE_EDM_FILTER_DEFINED_H

#include <complex>
#include "Parameters.hpp"
#include <iostream>

namespace ACE{

class EDM_Filter{
public:
  double threshold;
  double candidate_threshold;
  double increment_threshold;
  double sparse_prop_threshold;
  double threshold_in_round; double max_in_round;

  int max_coeffs;
  double mc_max; double mc_min; int mc_Nsteps;

  bool passes( const std::complex<double> &val) const;
  bool candidate_passes( const std::complex<double> &val) const;
  bool increment_passes( const std::complex<double> &val) ;
  bool sparse_prop_passes( const std::complex<double> &val) const;
  
  inline void new_round(){ 
    max_in_round=0.;
  }
 
  void print_info(std::ostream &ofs=std::cout)const;
  
  void setup(Parameters &param);

  EDM_Filter(double thr=0.) : threshold(thr){}
  EDM_Filter(Parameters &param){ 
    setup(param);
  }
};

}//namespace
#endif
