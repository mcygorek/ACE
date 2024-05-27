#ifndef ACE_BCF_DECOMPOSITION_MEIERTANNOR_DEFINED_H
#define ACE_BCF_DECOMPOSITION_MEIERTANNOR_DEFINED_H
#include "BCF_Decomposition.hpp"
#include "Parameters.hpp"

namespace ACE{

/*
J(omega) = (1/2) *  p*omega/[[(omega+Omega)^2+Gamma^2][(omega-Omega)^2+Gamma^2]]

see https://doi.org/10.1063/1.479669  
(There: J_there(omega) = pi * J_here(omega);
Note the symmetric limits in Eq. (8) )
*/

class BCF_Decomposition_MeierTannor: public BCF_Decomposition{
public:

  double p, Omega, Gamma, beta;

  //spectral density
  virtual double J(double omega)const;

  //non-Matsubara complex exponentials
  virtual std::vector<Term> get_nonMatsubara()const;
  //Matsubara complex exponentials
  virtual std::vector<Term> get_Matsubara(int N)const;

  virtual void print_info(std::ostream &ofs=std::cout)const;


  void setup(Parameters &param, const std::string &prefix="");
  BCF_Decomposition_MeierTannor(Parameters &param, const std::string &prefix=""){
    setup(param,prefix);
  }
  BCF_Decomposition_MeierTannor(){
    Parameters param;
    setup(param);
  }
};//class

}//namespace
#endif
