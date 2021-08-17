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

/*
  std::complex<double> integrate_times_expI(double xa, double xb, int Ndiscr, double arg) const{
  // int dx f(x) * exp(i*arg*x) 
  // single interval: int dx f(x) * exp(i*arg*x) -> (f(x2)+f(x1))/2 * dexp(i*arg*x)
  // dexp(i*arg*x) = int dx exp(i*arg*x) 
  //               = 1/(i*arg) *(exp(i*arg*x2) - exp(i*arg*x1))  
  //               = exp(i*arg*x1) * 1/(i*arg) (exp(i*arg*(x2-x1))-1)  
    
    if(fabs(arg)<1e-12)return integrate(xa, xb, Ndiscr);
  
    double h=(xb-xa)/(Ndiscr);
    std::complex<double> ep=exp(std::complex<double>(0.,arg*xa));
    std::complex<double> edp=exp(std::complex<double>(0.,arg*h));
    std::complex<double> dedp=std::complex<double>(0.,-1./arg)*
                              (exp(std::complex<double>(0.,arg*h))-1.);
 
    std::complex<double> res=0;
    T oldf2=f(xa); //is new f1
    for(int i=1; i<Ndiscr; i++){
      T newf2=f(xa+h*i);
      res+=(oldf2+newf2)/2.*ep*dedp;
      ep*=edp;
      oldf2=newf2;
    }
    return res;
  }
*/
  std::complex<double> integrate_times_expI(double xa, double xb, int Ndiscr, double arg) const{
    if(fabs(arg)<1e-12)return integrate(xa, xb, Ndiscr);
    int N=(Ndiscr+1)/2; //Nr. of intervals
    double h=(xb-xa)/(N);
    double h2=h/2.;
    double ah2=arg*h2;

    std::complex<double> ep = exp(std::complex<double>(0.,arg*(xa+h2)));
    std::complex<double> edp = exp(std::complex<double>(0.,arg*h));

    std::complex<double> A0 = (2./arg) * sin(ah2);
    std::complex<double> A1 = std::complex<double>(0.,1.) * (2./(arg*arg))
                              * (sin(ah2)-ah2*cos(ah2));
    std::complex<double> A2 = (2./(arg*arg*arg)) * 
                              ( (ah2*ah2-2.)*sin(ah2) + 2.*ah2*cos(ah2)) ;

    std::complex<double> res = 0.;
    T f1=f(xa);
    for(int n=0; n<N; n++){
      res += ep*f1*(A2-h2*A1)/(2.*h2*h2);
      res += ep*f(xa+(n+0.5)*h)*(A0-A2/(h2*h2));
        f1 = f(xa+(n+1.)*h);
      res += ep*f1*(A2+h2*A1)/(2.*h2*h2);
      ep *= edp;
    }
    return res;
  }


  virtual void print(const std::string &str, double xa, double xb, int Ndiscr){
    std::ofstream ofs(str.c_str());
    double h=(xb-xa)/Ndiscr;
    for(int i=0; i<=Ndiscr; i++){
      double x=xa+i*h;
      ofs<<x<<" "<<f(x)<<std::endl;
    }
  }
  virtual void print_integral(const std::string &str, double xa, double xb, int Ndiscr){
    std::ofstream ofs(str.c_str());
    double h=(xb-xa)/Ndiscr;
    T res=0;
    for(int i=0; i<Ndiscr; i++){
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
  for(int i=0; i<=Ndiscr; i++){
    double x=xa+i*h;
    std::complex<double> c=f(x);
    ofs<<x<<" "<<c.real()<<" "<<c.imag()<<std::endl;
  }
}
template<> void ComplexFunction::print_integral(const std::string &str, double xa, double xb, int Ndiscr){
  std::ofstream ofs(str.c_str());
  double h=(xb-xa)/Ndiscr;
  std::complex<double> res=0;
  for(int i=0; i<Ndiscr; i++){
    double x=xa+i*h;
    res+=h*0.5*(f(x)+f(x+h));
    ofs<<x<<" "<<res.real()<<" "<<res.imag()<<std::endl;
  }
}



typedef Smart_Ptr<RealFunction> RealFunctionPtr;
typedef Smart_Ptr<ComplexFunction> ComplexFunctionPtr;

//---- Some parameter independent instances of Function primitives
class RealFunction_Zero: public RealFunction{
public:
  virtual double f(double x)const{return 0.;}
};
RealFunctionPtr RealFunctionPtr_Zero=new RealFunction_Zero;

class RealFunction_Abs: public RealFunction{
public:
  virtual double f(double x)const{return fabs(x);}
  RealFunction_Abs(){}
};
RealFunctionPtr RealFunctionPtr_Abs=new RealFunction_Abs;



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

  virtual double f(double x)const{
    double y=(x-c)/width;
    if(y>=0.)return 1./(1.+exp(-y));
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

#endif
