#include "PCH.hpp"
#include <complex>
#include <fstream>
#include "Function.hpp"

namespace ACE{

template<typename T> 
  T Function_Interface<T>::integrate_Riemann(double xa, double xb, int Ndiscr)const{
    double h=(xb-xa)/(Ndiscr);
    T res=0.;
    for(int i=0; i<Ndiscr; i++){
      res+=f(xa+i*h)*h;
    }
    return res;
  }
template<typename T>
  T Function_Interface<T>::integrate(double xa, double xb, int Ndiscr)const{
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
template <typename T>
  std::complex<double> Function_Interface<T>::integrate_times_expI(double xa, double xb, int Ndiscr, double arg) const{
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

template <typename T>
  void Function_Interface<T>::print(const std::string &str, double xa, double xb, int Ndiscr){
    std::ofstream ofs(str.c_str());
    double h=(xb-xa)/Ndiscr;
    for(int i=0; i<=Ndiscr; i++){
      double x=xa+i*h;
      ofs<<x<<" "<<f(x)<<std::endl;
    }
  }
template<typename T>
  void Function_Interface<T>::print_integral(const std::string &str, double xa, double xb, int Ndiscr){
    std::ofstream ofs(str.c_str());
    double h=(xb-xa)/Ndiscr;
    T res=0;
    for(int i=0; i<Ndiscr; i++){
      double x=xa+i*h;
      res+=h*0.5*(f(x)+f(x+h));
      ofs<<x<<" "<<res<<std::endl;
    }
  }
 


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

template class Function_Interface<double>;
template class Function_Interface<std::complex<double> >;


//---- Some parameter independent instances of Function primitives
RealFunctionPtr RealFunctionPtr_Zero=std::make_shared<RealFunction_Zero>();

RealFunctionPtr RealFunctionPtr_Abs=std::make_shared<RealFunction_Abs>();



//---- Some parameter dependent standard functions


//-------- Compositions:

}//namespace
