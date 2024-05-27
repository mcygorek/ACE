#ifndef ACE_BCF_DECOMPOSITION_DRUDELORENTZ_DEFINED_H
#define ACE_BCF_DECOMPOSITION_DRUDELORENTZ_DEFINED_H
#include "BCF_Decomposition.hpp"
#include "Parameters.hpp"

namespace ACE{

/*
J(omega) = (2/pi) * lambda*gamma*omega/(gamma^2+omega^2)

see 10.1103/PhysRevResearch.5.013181 
*/

class BCF_Decomposition_DrudeLorentz: public BCF_Decomposition{
public:

  double lambda, gamma, beta;

  //spectral density
  virtual double J(double omega)const;

  //non-Matsubara complex exponentials
  virtual std::vector<Term> get_nonMatsubara()const;
  //Matsubara complex exponentials
  virtual std::vector<Term> get_Matsubara(int N)const;

  virtual void print_info(std::ostream &ofs=std::cout)const;

  //Terminator required
  virtual bool use_terminator()const;
  virtual std::complex<double> get_terminator(int Nk)const;

  void setup(Parameters &param, const std::string &prefix="");
  BCF_Decomposition_DrudeLorentz(Parameters &param, const std::string &prefix=""){
    setup(param,prefix);
  }
  BCF_Decomposition_DrudeLorentz(){
    Parameters param;
    setup(param);
  }
};//class

}//namespace
#endif
