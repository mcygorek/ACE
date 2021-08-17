#ifndef ACE_ENERGY_RANGE_DEFINED_H
#define ACE_ENERGY_RANGE_DEFINED_H

#include "Parameters.h"

//Define range of discretization
struct EnergyRange{
  std::pair<double, double> pair; 

  double E_min()const{ return pair.first; }
  double E_max()const{ return pair.second; }
  double E_range()const{ return pair.second-pair.first;}
  double E_mid()const{ return 0.5*(pair.second+pair.first);}
  double omega_min()const{ return pair.first/Constants::hbar_in_meV_ps; }
  double omega_max()const{ return pair.second/Constants::hbar_in_meV_ps; }
  
  std::string add_name(const std::string &mpgname, const std::string &str){
    if(mpgname=="")return str;
    return std::string(mpgname+"_"+str);
  }  
  void setup(Parameters &param, const std::string &mpgname){
    pair.first= Constants::hbar_in_meV_ps * 
                param.get_as_double(add_name(mpgname,"omega_min"));
    pair.first=param.get_as_double(add_name(mpgname,"E_min"),pair.first);

    pair.second= Constants::hbar_in_meV_ps * 
                 param.get_as_double(add_name(mpgname,"omega_max"));
    pair.second=param.get_as_double(add_name(mpgname,"E_max"),pair.second);

  }

  bool operator==(const EnergyRange &other){
    double error=1e-8;
    if(fabs(pair.first-other.pair.first)>error)return false;
    if(fabs(pair.second-other.pair.second)>error)return false;
    return true;
  }
  EnergyRange(Parameters &param, const std::string &mpgname=""){
    setup(param, mpgname);
  }
  EnergyRange(double Emax){ 
    pair.first=0; pair.second=Emax;
  }
  EnergyRange(double Emin, double Emax){ 
    pair.first=Emin; pair.second=Emax;
  }
  EnergyRange(){
    pair.first=0; pair.second=0;
  }
};


#endif
