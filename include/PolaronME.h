#ifndef POLARON_ME_DEFINED_H
#define POLARON_ME_DEFINED_H

#include "Function.hpp"
#include "Constants.hpp"

namespace ACE{

class PolaronME{
public:

  static std::complex<double> phi(double tau, RealFunctionPtr J, double temperature, double Emax=20., int Ndiscr=1e5){

    class f_cl_: public ComplexFunction{
      public: 
      RealFunctionPtr J;
      double T;
      double tau;
      virtual std::complex<double> f(double x) const{
        if(fabs(x)<1e-6)return 0.;
        double C=1.;
        if(T>1e-8){
          C=coth(hbar_in_meV_ps*x/(2.*kB_in_meV_by_K*T));
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
   
  //with shift_offset=(system eigenenergy difference):
  //   produces the _negative_ of the lamb shift. 
  static double polaron_shift(RealFunctionPtr J, double Emin=0., double Emax=20., int Ndiscr=1e5, double shift_offset=0.){
    
    class f_cl_: public RealFunction{
      public: 
      RealFunctionPtr J;
      double offset;
      virtual double f(double x) const{
        if(fabs(x-offset)<1e-6)return 0.;
        return J->f(x)/(x-offset);
      }
      f_cl_(RealFunctionPtr J_, double offs=0.): J(J_), offset(offs){}
    }f_cl(J, shift_offset);
  
    return hbar_in_meV_ps * f_cl.integrate(Emin/hbar_in_meV_ps, Emax/hbar_in_meV_ps, Ndiscr);
  }


};

}//namespace
#endif
