#ifndef ACE_SPECTRAL_DENSITY_SELECTOR_DEFINED_H
#define ACE_SPECTRAL_DENSITY_SELECTOR_DEFINED_H

#include <iosfwd>
#include "SpectralDensity.hpp"
#include "Function.hpp"

namespace ACE{

class SpectralDensity_Selector{
public:
  RealFunctionPtr SD;

  inline operator RealFunctionPtr() { return SD; } 

  static std::string add_name(const std::string & prefix, const std::string & str);
  
  static std::string get_prefix_without_J(const std::string & prefix="");
  
  static bool is_specified(Parameters &param, const std::string & prefix="");

  void setup(Parameters &param, const std::string & prefix="");
  
  SpectralDensity_Selector();
  
  inline SpectralDensity_Selector(Parameters &param, const std::string & prefix=""){
    setup(param, prefix);
  }
};

}//namespace
#endif
