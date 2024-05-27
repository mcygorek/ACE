#pragma once
#ifndef CONSTANTS_DEFINED_H
#define CONSTANTS_DEFINED_H
#include <cmath>
#include <complex>

namespace ACE{

constexpr double kB_in_eV_by_K=8.6173303e-5;       // [eV/K]
constexpr double kB_in_meV_by_K=8.6173303e-2;       // [meV/K]

constexpr double hbar_in_eV_ps=0.6582119569e-3;
constexpr double hbar_in_meV_ps=0.6582119569;
constexpr double inv_ps_to_meV=hbar_in_meV_ps; 
constexpr double meV_to_inv_ps=1./hbar_in_meV_ps; 

// inv_cm: c=299792458*100 cm/s -> hbar*omega=hbar*2*pi*nu=hbar*2*pi*c <-lambda=1cm
const double inv_cm_in_meV=0.12398409369863899086;
const double inv_cm_in_inv_ps=0.188364997625644;

const double muB_in_eV_by_T=5.7883818012e-5;   // [eV/T]
const double muB_in_meV_by_T=5.7883818012e-2;   // [meV/T]

const double eV_to_J=1.602176634e-19;
const double J_to_eV=1./eV_to_J;               // J=kg*m^2/s^2
const double J_to_meV=1000.*J_to_eV;

const double kg_to_meV_ps2_by_nm2=J_to_meV*1.0e6;      //kg=J*s^2/m^2

const double c_in_m_by_s=299792458; 
const double c_in_nm_by_ps=299792.458;

inline double beta_from_T(double T){
  if(T>1e12)return 0;
  if(T<1e-12)T=1e-12;
  return 1./(kB_in_meV_by_K*T);
}
inline double coth(double x){
  if(x<0.){
    return -coth(x);
  }
  if(fabs(x)<1e-12){
    return 0.; 
  }
  if(x<1e-6){
    double x2=x*x;
    return 1./x+x/3.-x2*x/45.+x2*x2*x*2./945.;
  }
  if(x>1e6){
    return -(1.+2./(exp(-2.*x)-1));
  }
  return 1.+2./(exp(2.*x)-1);
}
inline double gauss(double x, double s){
  double y=x/s;
  return exp(-0.5*y*y)/(s*sqrt(2.*M_PI));
}
inline double gauss_from_FWHM(double x, double fwhm){
  return gauss(x,fwhm/(2.*sqrt(2.*log(2.))));
}

inline double bose(double x){  // x=(E-mu)/(k_B T)
  return 1./(exp(x)-1.);
}
inline double fermi(double x){  // x=(E-mu)/(k_B T)
  if(x>1e4)return 0.;
  else if(x<-1e4)return 1.;
  else return 1./(exp(x)+1.);
}

inline int factorial(int n, int stop=1){
  if(n<=stop)return 1;
  return factorial(n-1)*n;
}

inline double logistic(double x){
  if(x>=0.)return 1./(1.+exp(-x));
  else return 1.-1./(1.+exp(x));
}

inline std::complex<double> complex_coth(std::complex<double> arg){
  if(abs(arg)<1e-6){
    return 1./arg+arg/3.-arg*arg*arg/45.;
  }else if(real(arg)>1.){
    return (1.+exp(-2.*arg))/(1.-exp(-2.*arg));
  }else{
    return (exp(2.*arg)+1.)/(exp(2.*arg)-1.);
  }
}

inline std::complex<double> complex_cot(std::complex<double> arg){
  return 1./tan(arg);
}
};


#endif
