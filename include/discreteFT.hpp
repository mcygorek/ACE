#pragma once
#ifndef DISCRETEFT_DEFINED_H
#define DISCRETEFT_DEFINED_H

#include <complex> 
#include <vector> 
#include "Simulation_Results.hpp"

namespace ACE{


//DFT of {in[0], in[stride], in[2*stride], ..., in[(N-1)*stride }
void discreteFT(std::complex<double> *in, std::complex<double> *out, int N, int stride, int sign=-1);

//zero-pads input to next power of 2
std::vector<std::complex<double> > discreteFT(
     const std::vector<std::complex<double> > & in, int sign=-1);

//Fourier integral a column of Simulation_Results
Simulation_Results resultsFFT(const Simulation_Results & in, int col, int sign=-1, int integrate_mode=0);

//correction terms for Fourier integrals (cf. "Numerical Recipes" 13.9)
double FFT_trapezoidal_correction_W(double t);
std::complex<double> FFT_trapezoidal_correction_a0(double);


} //namespace
#endif
