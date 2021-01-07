#ifndef FUNCTION_DEFINED_H_
#define FUNCTION_DEFINED_H_

#include <complex>
#include <fstream>
#include "Smart_Ptr.h"

template <typename T> class Function_Interface{
public:

  virtual T f(double x) const=0;
 
  T integrate_Riemann(double xa, double xb, int Ndiscr) const{
    double h=(xb-xa)/(Ndiscr);
    T res=0.;
    for(int i=0; i<Ndiscr; i++){
      res+=f(xa+i*h)*h;
    } 
    return res;
  } 
  ///integrate using Simpson's rule
  T integrate(double xa, double xb, int Ndiscr) const{
    int N=(Ndiscr+1)/2;
    double h=(xb-xa)/(N);
    double h6=h/6.;

    T res=h6*(f(xa)+f(xb));
    for(int i=1; i<N; i++){
      res+=h6*2.*f(xa+i*h);
    } 
    for(int i=1; i<=N; i++){
      res+=h6*4.*f(xa+(i-0.5)*h);
    } 
    return res;
  }

  void print(const std::string &str, double xa, double xb, int Ndiscr){
    std::ofstream ofs(str.c_str());
    double h=(xb-xa)/Ndiscr;
    for(size_t i=0; i<=Ndiscr; i++){
      double x=xa+i*h;
      ofs<<x<<" "<<f(x)<<std::endl;
    }
  }
  void print_integral(const std::string &str, double xa, double xb, int Ndiscr){
    std::ofstream ofs(str.c_str());
    double h=(xb-xa)/Ndiscr;
    T res=0;
    for(size_t i=0; i<Ndiscr; i++){
      double x=xa+i*h;
      res+=h*0.5*(f(x)+f(x+h));
      ofs<<x<<" "<<res<<std::endl;
    }
  }
 
  virtual ~Function_Interface(){}
};

typedef Function_Interface<double> RealFunction;
typedef Function_Interface<std::complex<double> > ComplexFunction;

template<> void ComplexFunction::print(const std::string &str, double xa, double xb, int Ndiscr){
  std::ofstream ofs(str.c_str());
  double h=(xb-xa)/Ndiscr;
  for(size_t i=0; i<=Ndiscr; i++){
    double x=xa+i*h;
    std::complex<double> c=f(x);
    ofs<<x<<" "<<c.real()<<" "<<c.imag()<<std::endl;
  }
}
template<> void ComplexFunction::print_integral(const std::string &str, double xa, double xb, int Ndiscr){
  std::ofstream ofs(str.c_str());
  double h=(xb-xa)/Ndiscr;
  std::complex<double> res=0;
  for(size_t i=0; i<Ndiscr; i++){
    double x=xa+i*h;
    res+=h*0.5*(f(x)+f(x+h));
    ofs<<x<<" "<<res.real()<<" "<<res.imag()<<std::endl;
  }
}



typedef Smart_Ptr<RealFunction> RealFunctionPtr;
typedef Smart_Ptr<ComplexFunction> ComplexFunctionPtr;


class RealFunction_Zero_Class: public RealFunction{
  double f(double x)const{return 0.;}
}RealFunction_Zero;

RealFunctionPtr RealFunctionPtr_Zero=new RealFunction_Zero_Class;

#endif
