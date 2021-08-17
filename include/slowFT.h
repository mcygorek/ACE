#ifndef SLOWFT_DEFINED_H
#define SLOWFT_DEFINED_H

/** Slow Fourier transform of an interpolated function. For test purposes. */

#include <vector>
#include <complex>
#include <cmath>
#include <fstream>
#include "Constants.h"


std::vector<std::complex<double> > slowFT_Simpson(const std::vector<std::complex<double> > & in, double wa, double wb, int Ndiscr, double ta, double dt, int sig=1, int Nsubdiv=10){

  if(Nsubdiv<1){
    std::cerr<<"Error: slowFT: Nsubdiv<1!"<<std::endl;
    exit(1);
  }

  std::vector<std::complex<double> > ret(Ndiscr+1);
  double dw=(wb-wa)/Ndiscr;
  double dt2=dt/(2.*Nsubdiv);

  for(int i=0; i<=Ndiscr; i++){
    double w=wa+i*dw;
#ifdef SLOWFT_PRINT
    std::cout<<"slowFT: E="<<w*Constants::hbar_in_meV_ps<<std::endl;
#endif 
    std::complex<double> res=0;
    for(int j=0; j<(int)in.size()-1; j++){
      double t0=ta+j*dt;
      std::complex<double> dy=(in[j+1]-in[j])/(2.*Nsubdiv);

      std::complex<double> eph=exp(std::complex<double>(0, sig*w*t0));
      std::complex<double> deph=exp(std::complex<double>(0, sig*w*dt2));

      res+=1./3.*dt2*in[j]*eph;
//      res+=1./3.*dt2*in[j]*exp(std::complex<double>(0, sig*w*t0));

      for(int k=0; k<Nsubdiv; k++){
        double ddt=(2.*k+1.)*dt2;
//        double t=t0+ddt;
        std::complex<double> y=(in[j]+dy*ddt);
        eph*=deph;
        res+= 4./3.*dt2 * eph * y;
//        res+= 4./3.*dt2 * exp(std::complex<double>(0, sig*w*t)) * y;
        if(k<Nsubdiv-1){
          eph*=deph;
          res+= 2./3.*dt2 * eph * (y+dy*dt2);
//          res+= 2./3.*dt2 * exp(std::complex<double>(0, sig*w*(t+dt2))) * (y+dy*dt2);
        }
      }
      eph*=deph;
      res+=1./3.*dt2*in[j+1]*eph;
//      res+=1./3.*dt2*in[j+1]*exp(std::complex<double>(0, sig*w*(t0+dt)));
    }
    ret[i]=res;
  }

  return ret;
}
std::vector<std::complex<double> > slowFT_Riemann(const std::vector<std::complex<double> > & in, double wa, double wb, int Ndiscr, double ta, double dt, int sig=1, int Nsubdiv=10){


  std::vector<std::complex<double> > ret(Ndiscr+1);
  double dw=(wb-wa)/Ndiscr;
  double dt2=dt/Nsubdiv;

  for(int i=0; i<=Ndiscr; i++){
    double w=wa+i*dw;
    std::complex<double> res=0;
    for(int j=0; j<(int)in.size()-1; j++){
      for(int k=0; k<Nsubdiv; k++){
        double t=ta+j*dt+k*dt2;
        std::complex<double> y=((double)(Nsubdiv-k)*in[j]+(double)k*in[j+1])/((double)(Nsubdiv));
        res+= dt2 * exp(std::complex<double>(0, sig*w*t)) * y;
      }
    }
    ret[i]=res;
  }

  return ret;
}

std::vector<std::complex<double> > slowFT(const std::vector<std::complex<double> > & in, double wa, double wb, int Ndiscr, double ta, double dt, int sig=1, int Nsubdiv=10){
#ifdef SLOWFT_USE_RIEMANN
  return slowFT_Riemann(in, wa, wb, Ndiscr, ta, dt, sig, Nsubdiv);
#else 
  return slowFT_Simpson(in, wa, wb, Ndiscr, ta, dt, sig, Nsubdiv);
#endif
}

void print_slowFT(const std::string &fname, const std::vector<std::complex<double> > & in, double wa, double wb, int Ndiscr, double ta, double dt, int sig=1, int Nsubdiv=10){

  std::vector<std::complex<double> > out=slowFT(in, wa, wb, Ndiscr, ta, dt, sig, Nsubdiv);

  std::ofstream ofs(fname.c_str());
  for(int i=0; i<(int)out.size(); i++){
    ofs<<wa+i*(wb-wa)/Ndiscr<<" "<<out[i].real()<<" "<<out[i].imag()<<std::endl;
  }

}
void print_slowFT_meV(const std::string &fname, const std::vector<std::complex<double> > & in, double Ea, double Eb, int Ndiscr, double ta, double dt, int sig=1, int Nsubdiv=10){

  double hbar=Constants::hbar_in_meV_ps;

  std::vector<std::complex<double> > out=slowFT(in, Ea/hbar, Eb/hbar, Ndiscr, ta, dt, sig, Nsubdiv);

  std::ofstream ofs(fname.c_str());
  for(int i=0; i<(int)out.size(); i++){
    ofs<<Ea+i*(Eb-Ea)/Ndiscr<<" "<<out[i].real()<<" "<<out[i].imag()<<std::endl;
  }

}
#endif
