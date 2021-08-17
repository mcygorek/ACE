#ifndef CONSTANTS_DEFINED_H
#define CONSTANTS_DEFINED_H
#include <cmath>

namespace Constants{

const double kB_in_eV_by_K=8.6173303e-5;       // [eV/K]
const double kB_in_meV_by_K=8.6173303e-2;       // [meV/K]

const double hbar_in_eV_ps=0.6582119569e-3;
const double hbar_in_meV_ps=0.6582119569;
const double inv_ps_to_meV=hbar_in_meV_ps; 
const double meV_to_inv_ps=1./hbar_in_meV_ps; 

const double muB_in_eV_by_T=5.7883818012e-5;   // [eV/T]
const double muB_in_meV_by_T=5.7883818012e-2;   // [meV/T]

const double eV_to_J=1.602176634e-19;
const double J_to_eV=1./eV_to_J;               // J=kg*m^2/s^2
const double J_to_meV=1000.*J_to_eV;

const double kg_to_meV_ps2_by_nm2=J_to_meV*1.0e6;      //kg=J*s^2/m^2

const double c_in_m_by_s=299792458; 
const double c_in_nm_by_ps=299792.458;

double beta_from_T(double T){
  if(T>1e12)return 0;
  if(T<1e-12)T=1e-12;
  return 1./(kB_in_meV_by_K*T);
}
double coth(double x){
  if(fabs(x)<1e-10) return 0.; 
  if(x<0.)return -coth(x);
  if(x>1e6) return 1.;
  return 1.+2./(exp(2.*x)-1);
}
double gauss(double x, double s){
  double y=x/s;
  return exp(-0.5*y*y)/(s*sqrt(2.*M_PI));
}
double gauss_from_FWHM(double x, double fwhm){
  return gauss(x,fwhm/(2.*sqrt(2.*log(2.))));
}

double bose(double x){  // x=(E-mu)/(k_B T)
  return 1./(exp(x)-1.);
}
double fermi(double x){  // x=(E-mu)/(k_B T)
  if(x>1e4)return 0.;
  else if(x<-1e4)return 1.;
  else return 1./(exp(x)+1.);
}

int factorial(int n, int stop=1){
  if(n<=stop)return 1;
  return factorial(n-1)*n;
}

double logistic(double x){
  if(x>=0.)return 1./(1.+exp(-x));
  else return 1.-1./(1.+exp(x));
}

};


#endif
