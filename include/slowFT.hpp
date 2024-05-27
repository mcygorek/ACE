#pragma once
#ifndef SLOWFT_DEFINED_H
#define SLOWFT_DEFINED_H

/** Slow Fourier transform of an interpolated function. For test purposes. */

#include <vector>
#include <complex>

namespace ACE{

extern std::vector<std::complex<double> > slowFT_Simpson(const std::vector<std::complex<double> > & in, double wa, double wb, int Ndiscr, double ta, double dt, int sig=1, int Nsubdiv=10);

extern std::vector<std::complex<double> > slowFT_Riemann(const std::vector<std::complex<double> > & in, double wa, double wb, int Ndiscr, double ta, double dt, int sig=1, int Nsubdiv=10);

inline std::vector<std::complex<double> > slowFT(const std::vector<std::complex<double> > & in, double wa, double wb, int Ndiscr, double ta, double dt, int sig=1, int Nsubdiv=10){
#ifdef SLOWFT_USE_RIEMANN
  return slowFT_Riemann(in, wa, wb, Ndiscr, ta, dt, sig, Nsubdiv);
#else 
  return slowFT_Simpson(in, wa, wb, Ndiscr, ta, dt, sig, Nsubdiv);
#endif
}

extern void print_slowFT(const std::string &fname, const std::vector<std::complex<double> > & in, double wa, double wb, int Ndiscr, double ta, double dt, int sig=1, int Nsubdiv=10);

extern void print_slowFT_meV(const std::string &fname, const std::vector<std::complex<double> > & in, double Ea, double Eb, int Ndiscr, double ta, double dt, int sig=1, int Nsubdiv=10);


}//namespace
#endif
