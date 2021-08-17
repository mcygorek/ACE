#ifndef POLARON_ME_DEFINED_H
#define POLARON_ME_DEFINED_H


#include "Function.h"
#include "Constants.h"

class PolaronME{
public:

  static std::complex<double> phi(double tau, RealFunctionPtr J, double temperature, double Emax=20., int Ndiscr=1e5){

    class f_cl_: public ComplexFunction{
      public: 
      RealFunctionPtr J;
      double T;
      double tau;
      virtual std::complex<double> f(double x) const{
        if(x<1e-6)return 0.;
        double C=1.;
        if(T>1e-8){
          C=Constants::coth(Constants::hbar_in_meV_ps*x)/(2.*Constants::kB_in_meV_by_K*T);
        }
        return J->f(x)/(x*x)* std::complex<double>(C*cos(x*tau), -sin(x*tau));
      }
      f_cl_(double tau_, RealFunctionPtr J_, double T_): J(J_), T(T_), tau(tau_){}
    }f_cl(tau, J, temperature);
  
    return f_cl.integrate(0, Emax, Ndiscr);
  }


  static double renorm(RealFunctionPtr J, double temperature, double Emax=20., int Ndiscr=1e5){
    return exp(-phi(0., J, temperature, Emax, Ndiscr).real()/2.);
  } 
   

  static double polaron_shift(RealFunctionPtr J, double Emax=20., int Ndiscr=1e5){

    class f_cl_: public RealFunction{
      public: 
      RealFunctionPtr J;
      virtual double f(double x) const{
        if(x<1e-6)return 0.;
        return J->f(x)/x;
      }
      f_cl_(RealFunctionPtr J_): J(J_){}
    }f_cl(J);
  
    return Constants::hbar_in_meV_ps * f_cl.integrate(0, Emax/Constants::hbar_in_meV_ps, Ndiscr);
  }


};

#endif
