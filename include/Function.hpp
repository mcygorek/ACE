#pragma once
#ifndef FUNCTION_DEFINED_H_
#define FUNCTION_DEFINED_H_

#include <complex>
#include <memory>
#include <fstream>

namespace ACE{

template <typename T> class Function_Interface{
public:

  virtual T f(double x) const=0;
 
  T integrate_Riemann(double xa, double xb, int Ndiscr)const;
   
  ///integrate using Simpson's rule
  T integrate(double xa, double xb, int Ndiscr)const;

  std::complex<double> integrate_times_expI(double xa, double xb, int Ndiscr, double arg) const;

  virtual void print(const std::string &str, double xa, double xb, int Ndiscr);
  
  virtual void print_integral(const std::string &str, double xa, double xb, int Ndiscr);

  virtual ~Function_Interface(){}
};

typedef Function_Interface<double> RealFunction;
typedef Function_Interface<std::complex<double> > ComplexFunction;

template<> void ComplexFunction::print(const std::string &str, double xa, double xb, int Ndiscr);

template<> void ComplexFunction::print_integral(const std::string &str, double xa, double xb, int Ndiscr);


typedef std::shared_ptr<RealFunction> RealFunctionPtr;
typedef std::shared_ptr<ComplexFunction> ComplexFunctionPtr;


//---- Some parameter independent instances of Function primitives
class RealFunction_Zero: public RealFunction{
public:
  virtual double f(double x)const{return 0.;}
};
extern RealFunctionPtr RealFunctionPtr_Zero;

class RealFunction_Abs: public RealFunction{
public:
  virtual double f(double x)const{return fabs(x);}
  RealFunction_Abs(){}
};
extern RealFunctionPtr RealFunctionPtr_Abs;



//---- Some parameter dependent standard functions

class RealFunction_Const: public RealFunction{
public:
  double a;
  virtual double f(double x)const{return a;}
  RealFunction_Const(double a_=0.) : a(a_){}
};


class RealFunction_Lorentzian: public RealFunction{
public:
  double scale, gamma, c;

  virtual double f(double x)const{
    return scale/M_PI*gamma/( (x-c)*(x-c) + gamma*gamma );
  }
  RealFunction_Lorentzian(double scale_=1., double gamma_=0.01, double c_=0.){
    scale=scale_;
    gamma=gamma_;
    c=c_;
  }
};

class RealFunction_Logistic: public RealFunction{
public:
  double c, width;
  // f(x)= 1/(1+exp(-(x-c)/w)
 
  virtual double f(double x)const{
    if(fabs(width)<1e-16){
      if(x==c)return 0.5;
      if(width>=0.){
        return (x>c) ? 1. : 0.;
      }else{
        return (x>c) ? 0. : 1.;
      }
    }
    double y=(x-c)/width;
    if(fabs(y)<1e-5)return 0.5;
    if(y>0.) return 1./(1.+exp(-y));
    else return 1.-1./(1.+exp(y));
  }
  RealFunction_Logistic(double c_=0., double width_=1.)
   : c(c_), width(width_){
  }
};


//-------- Compositions:

class RealFunction_Scale: public RealFunction{
public:
  RealFunctionPtr p;
  double scale;
  virtual double f(double x)const{
    return scale*p->f(x);
  }

  RealFunction_Scale(){}
  RealFunction_Scale(RealFunctionPtr &p_, double scale_) 
   : p(p_), scale(scale_){
  }
};
//Note the minus sign in shift: shifts the function, not the argument!
class RealFunction_ScaleShift: public RealFunction{
public:
  RealFunctionPtr p;
  double scale;
  double shift;

  virtual double f(double x)const{
    return scale*p->f(x-shift);
  }

  RealFunction_ScaleShift(){}
  RealFunction_ScaleShift(RealFunctionPtr &p_, double scale_, double shift_) 
   : p(p_), scale(scale_), shift(shift_){
  }
};

class RealFunction_Product: public RealFunction{
public:
  RealFunctionPtr F1, F2;
  virtual double f(double x)const{
    return F1->f(x) * F2->f(x);
  }
  RealFunction_Product(RealFunctionPtr &F1_, RealFunctionPtr &F2_)
   : F1(F1_), F2(F2_){
  }
};

class RealFunction_Sum: public RealFunction{
public:
  RealFunctionPtr F1, F2;
  virtual double f(double x)const{
    return F1->f(x) + F2->f(x);
  }
  RealFunction_Sum(RealFunctionPtr &F1_, RealFunctionPtr &F2_)
   : F1(F1_), F2(F2_){
  }
};

class RealFunction_Chain: public RealFunction{
public:
  RealFunctionPtr F1, F2;
  virtual double f(double x)const{
    return F1->f( F2->f(x) );
  }
  RealFunction_Chain(RealFunctionPtr &F1_, RealFunctionPtr &F2_)
   : F1(F1_), F2(F2_){
  }
};

//extern template class Function_Interface<double>;
//extern template class Function_Interface<std::complex<double> >;

}//namespace
#endif
