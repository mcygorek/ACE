#ifndef SPECTRAL_DENSITY_DEFINED_H
#define SPECTRAL_DENSITY_DEFINED_H

#include "Function.hpp"
//#include <cstdlib>
//#include <iostream>
//#include "Constants.hpp"
//#include "Parameters.hpp"

namespace ACE{
class Parameters;

/** Spectral density for quantum dot  */
class SpectralDensity_QD: public RealFunction{
public:
  double mass_density;
  double a_e, a_h;
  double c_s;
  double D_e, D_h;

  std::string add_name(const std::string &mpgname, const std::string &str)const;
 
  void setup(Parameters &param, const std::string &prefix="");

  virtual double f(double x)const;

  SpectralDensity_QD();
  SpectralDensity_QD(Parameters &param, const std::string &prefix="");
  
  inline virtual ~SpectralDensity_QD(){}
};


/** (sub-/super-)ohmic spectral density  */
class SpectralDensity_sohmic: public RealFunction{
public:
  double s, alpha, omega_c;
  bool div_by_cutoff;
  enum CMODE { NONE, EXP, GAUSS, DRUDE, HARD } cutoff_mode;

  std::string add_name(const std::string &mpgname, const std::string &str)const;
 
  void setup(Parameters &param, const std::string &prefix="");

  virtual double f(double x)const;
  
  SpectralDensity_sohmic();
  
  SpectralDensity_sohmic(Parameters &param, const std::string &prefix="");

  virtual ~SpectralDensity_sohmic(){}
};

/** Bump function: smooth function with compact support */
class SpectralDensity_bump: public RealFunction{
public:
  double wmin,wmax,scale;

  std::string add_name(const std::string &mpgname, const std::string &str)const;

  void setup(Parameters &param, const std::string &prefix="");

  virtual double f(double x)const;
  
  SpectralDensity_bump(double wmin_=-1, double wmax_=1., double scale_=1.)
    : wmin(wmin_), wmax(wmax_), scale(scale_) {
  }
  
  SpectralDensity_bump(Parameters &param, const std::string &prefix="");

  virtual ~SpectralDensity_bump(){}
};

}//namespace
#endif
