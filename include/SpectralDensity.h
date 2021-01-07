#ifndef SPECTRAL_DENSITY_DEFINED_H
#define SPECTRAL_DENSITY_DEFINED_H

#include "Function.h"
#include <cstdlib>
#include <iostream>
#include "Constants.h"


/** Spectral density for quantum dot  */
class SpectralDensity_QD: public RealFunction{
public:
  double mass_density;
  double a_e, a_h;
  double c_s;
  double D_e, D_h;

  void set_defaults(){
    D_e=7.0*1000;
    D_h=-3.5*1000;
    mass_density=5370*Constants::kg_to_meV_ps2_by_nm2*1.0e-27;
    c_s=5110*1e-3;
    a_e=4;
    a_h=a_e/1.15;
  }

  virtual double f(double x)const{
    double x2=x*x;
    double c2=c_s*c_s;
    double arg=D_e*exp(-x2*a_e*a_e/(4.*c2))-D_h*exp(-x2*a_h*a_h/(4.*c2));
    return x2*x/(4.*M_PI*M_PI*mass_density*Constants::hbar_in_meV_ps*c2*c2*c_s)*arg*arg;
  }

  SpectralDensity_QD(){
    set_defaults();
  }

};



#endif





