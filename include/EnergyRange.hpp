#ifndef ACE_ENERGY_RANGE_DEFINED_H
#define ACE_ENERGY_RANGE_DEFINED_H

#include <utility>
#include <string>
//#include "Parameters.hpp"
#include "Constants.hpp"

namespace ACE{
class Parameters;

//Define range of discretization
struct EnergyRange{
  std::pair<double, double> pair; 

  inline double E_min()const{ return pair.first; }
  inline double E_max()const{ return pair.second; }
  inline double E_range()const{ return pair.second-pair.first;}
  inline double E_mid()const{ return 0.5*(pair.second+pair.first);}
  inline double omega_min()const{ return pair.first/hbar_in_meV_ps; }
  inline double omega_max()const{ return pair.second/hbar_in_meV_ps; }
  
  std::string add_name(const std::string &mpgname, const std::string &str)const;
    
  void setup(Parameters &param, const std::string &mpgname);

  bool operator==(const EnergyRange &other)const;
  
  inline EnergyRange(Parameters &param, const std::string &mpgname=""){
    setup(param, mpgname);
  }
  inline EnergyRange(double Emax){ 
    pair.first=0; pair.second=Emax;
  }
  inline EnergyRange(double Emin, double Emax){ 
    pair.first=Emin; pair.second=Emax;
  }
  inline EnergyRange(){
    pair.first=0; pair.second=0;
  }
};

}//namespace
#endif
