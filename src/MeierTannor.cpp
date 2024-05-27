#include "MeierTannor.hpp"
#include "Constants.hpp"

namespace ACE{

double MeierTannor::J(double omega)const{
  return eta* omega/(((omega+Omega)*(omega+Omega)+Gamma*Gamma)*((omega-Omega)*(omega-Omega)+Gamma*Gamma));
}
std::complex<double> MeierTannor::J(std::complex<double> omega)const{
  return eta* omega/(((omega+Omega)*(omega+Omega)+Gamma*Gamma)*((omega-Omega)*(omega-Omega)+Gamma*Gamma));
}

std::complex<double> MeierTannor::alpha_plus(double temperature)const{
  std::complex<double> coth=1.;
  if(temperature>1e-8){
    double beta=hbar_in_meV_ps/(kB_in_meV_by_K*temperature);
    std::complex<double> arg=(beta/2.)*std::complex<double>(Omega, Gamma);
    if(abs(arg)<1e-6){
      coth=1./arg+arg/3.-arg*arg*arg/45.;
    }else if(real(arg)>1.){
      coth=(1.+exp(-2.*arg))/(1.-exp(-2.*arg));
    }else{
      coth=(exp(2.*arg)+1.)/(exp(2.*arg)-1.);
    }
  }
  return 4.*eta/(Omega*Gamma)*(coth + std::complex<double>(0, 1.));
}
std::complex<double> MeierTannor::alpha_minus(double temperature)const{
  std::complex<double> coth=1.;
  if(temperature>1e-8){
    double beta=hbar_in_meV_ps/(kB_in_meV_by_K*temperature);
    std::complex<double> arg=(beta/2.)*std::complex<double>(Omega, -Gamma);
    if(abs(arg)<1e-6){
      coth=1./arg+arg/3.-arg*arg*arg/45.;
    }else if(real(arg)>1.){
      coth=(1.+exp(-2.*arg))/(1.-exp(-2.*arg));
    }else{
      coth=(exp(2.*arg)+1.)/(exp(2.*arg)-1.);
    }
  }
  return 4.*eta/(Omega*Gamma)*(coth - std::complex<double>(0, 1.));
}
std::complex<double> MeierTannor::alpha_therm(int k, double temperature)const{
  double beta=hbar_in_meV_ps/(kB_in_meV_by_K*temperature);
  double nu=2.*M_PI/beta;
  return std::complex<double>(0,2.)*nu*J(std::complex<double>(0,k*nu));
}


}
