#ifndef ACE_BCF_DECOMPOSITION_DEFINED_H
#define ACE_BCF_DECOMPOSITION_DEFINED_H
#include <utility>
#include <vector>
#include <complex>
#include <iostream>
	

namespace ACE{

/** A central part of HEOM approaches is to fit a given bath correlation 
function to a sum of (complex) exponentials. In general, this leads to an 
infinite series with a finite number of terms originating from poles of the
spectral density and an infinite number of Matsubara terms from poles of the
Bose (and/or Fermi?) distribution at finite temperatures.
(More terms are needed for smaller temperatures; T->0 ill-defined).
Refs:
[1]: 
[2]: Meier and Tannor: 10.1063/1.479669
[3]: QuTiP-BoFiN: 10.1103/PhysRevResearch.5.013181 


Our strategy is to calculate PTs for each term (single damped oscillator) and
sum them up ACE-style.

This class is an interface, from which one obtains the decomposition.
To be overloaded by classes implementing concrete spectral densities.
*/

class BCF_Decomposition{
public:

  //Coefficients: c*exp(-nu*t)  ->  (c, nu)
  typedef std::pair<std::complex<double>, std::complex<double> >  Term;

  //just for reference: spectral density
  virtual double J(double omega)const;

  //non-Matsubara complex exponentials
  virtual std::vector<Term> get_nonMatsubara()const = 0;
  //Matsubara complex exponentials
  virtual std::vector<Term> get_Matsubara(int N)const = 0;
  //first non-Matsubara, then Matsubara terms:
  virtual std::vector<Term> get_all(int N_matsubara)const;

  virtual void print_all(int N_matsubara, std::ostream &ofs=std::cout)const;
  virtual void print_info(std::ostream &ofs=std::cout)const=0;
  virtual void print_BCF(const std::string &fname, double dt, double te, int N_matsubara)const;

  virtual void print_J(const std::string &str, double xa, double xb, int Ndiscr)const;

  //Terminator required, e.g., for Drude-Lorentz (see [3])
  virtual bool use_terminator()const;
  virtual std::complex<double> get_terminator(int Nk)const;

};//class

}//namespace
#endif
