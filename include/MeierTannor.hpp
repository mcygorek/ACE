#ifndef ACE_MEIERTANNOR_DEFINED_H
#define ACE_MEIERTANNOR_DEFINED_H

#include <complex>
#include <cmath>

namespace ACE{

/** Spectral density of form: 
J(omega)=eta * omega/(((omega+Omega)^2+Gamma^2)((omega+Omega)^2+Gamma^2))
*/
class MeierTannor{
public: 
  double eta, Omega, Gamma;

  double J(double omega)const;
  std::complex<double> J(std::complex<double> omega)const;
 
  std::complex<double> alpha_plus(double temperature)const;
  std::complex<double> alpha_minus(double temperature)const;
  std::complex<double> alpha_therm(int k, double temperature)const;

  MeierTannor(double eta_, double Omega_, double Gamma_) 
    : eta(eta_), Omega(Omega_), Gamma(Gamma_)  {}
  MeierTannor(){}
};//class
}//namespace
#endif 
