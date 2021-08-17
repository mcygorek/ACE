#ifndef ACE_READ_TEMPERATURE_DEFINED_H_
#define ACE_READ_TEMPERATURE_DEFINED_H_

#include "Parameters.h"
#include "Constants.h"

double readTemperature(Parameters &param, const std::string &prefix=""){
  double temperature=0.;
  double fac=Constants::hbar_in_meV_ps/Constants::kB_in_meV_by_K;
  temperature=param.get_as_double("temperature",temperature);
  if(param.is_specified("temperature_unitless"))temperature=param.get_as_double("temperature_unitless",temperature/fac)*fac;

  if(prefix!=""){
    temperature=param.get_as_double(std::string(prefix+"_temperature"),temperature);
    if(param.is_specified(std::string(prefix+"_temperature_unitless")))temperature=param.get_as_double(std::string(prefix+"_temperature_unitless"),temperature/fac)*fac;
  }
 
  return temperature;
}


#endif
