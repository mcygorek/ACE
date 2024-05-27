#include "PCH.hpp"
#include "Parameters.hpp"
#include "Constants.hpp"
#include "ReadTemperature.hpp"
#include "ReaderBasics.hpp"

namespace ACE{

double readTemperature(Parameters &param, const std::string &prefix){
  double temperature=0.;
  double fac=hbar_in_meV_ps/kB_in_meV_by_K;
  temperature=param.get_as_double("temperature",temperature);
  if(param.is_specified("temperature_unitless"))temperature=param.get_as_double("temperature_unitless",temperature/fac)*fac;

  if(prefix!=""){
    temperature=param.get_as_double(std::string(prefix+"_temperature"),temperature);
    if(param.is_specified(std::string(prefix+"_temperature_unitless")))temperature=param.get_as_double(std::string(prefix+"_temperature_unitless"),temperature/fac)*fac;
  }
 
  return temperature;
}

}//namespace
