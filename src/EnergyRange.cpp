#include "EnergyRange.hpp"
#include "Parameters.hpp"
#include "Constants.hpp"

namespace ACE{
//Define range of discretization
  
  std::string EnergyRange::add_name(const std::string &mpgname, const std::string &str)const{
    if(mpgname=="")return str;
    return std::string(mpgname+"_"+str);
  }
  
  void EnergyRange::setup(Parameters &param, const std::string &mpgname){
    pair.first= hbar_in_meV_ps * 
                param.get_as_double(add_name(mpgname,"omega_min"));
    pair.first=param.get_as_double(add_name(mpgname,"E_min"),pair.first);

    pair.second= hbar_in_meV_ps * 
                 param.get_as_double(add_name(mpgname,"omega_max"));
    pair.second=param.get_as_double(add_name(mpgname,"E_max"),pair.second);

  }

  bool EnergyRange::operator==(const EnergyRange &other)const{
    double errorbar=1e-8;
    if(fabs(pair.first-other.pair.first)>errorbar)return false;
    if(fabs(pair.second-other.pair.second)>errorbar)return false;
    return true;
  }

}//namespace
