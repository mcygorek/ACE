#include "BCF_Decomposition_MeierTannor.hpp"
#include "ReadTemperature.hpp"
#include "Constants.hpp"
#include <fstream>


namespace ACE{

double BCF_Decomposition_MeierTannor::J(double omega)const{
  return 0.5*p*omega/(((omega+Omega)*(omega+Omega)+Gamma*Gamma)*((omega-Omega)*(omega-Omega)+Gamma*Gamma));
}


std::vector<BCF_Decomposition::Term> BCF_Decomposition_MeierTannor::get_nonMatsubara()const{
  std::vector<BCF_Decomposition::Term> terms(2);
  
//  terms[0].first=p/(Omega*Gamma)*(complex_coth(beta/2.*std::complex<double>(Omega, Gamma))+1.);
//  terms[1].first=p/(Omega*Gamma)*(complex_coth(beta/2.*std::complex<double>(Omega, -Gamma))-1.);
  terms[0].first=M_PI/16.*p/(Omega*Gamma)*(complex_coth(beta/2.*std::complex<double>(Omega, Gamma))-1.);
  terms[1].first=M_PI/16.*p/(Omega*Gamma)*(complex_coth(beta/2.*std::complex<double>(Omega, -Gamma))+1.);
  terms[0].second=std::complex<double>(Gamma,-Omega);
  terms[1].second=std::complex<double>(Gamma,Omega);
  
  return terms;
}
std::vector<BCF_Decomposition::Term> BCF_Decomposition_MeierTannor::get_Matsubara(int N)const{
  std::vector<BCF_Decomposition::Term> terms(N);
  for(int k=0; k<N; k++){
    double nu=2.*M_PI*(k+1)/beta;
    std::complex<double> denom(Omega*Omega-nu*nu+Gamma*Gamma, 2.*nu*Omega);
    terms[k].first=-M_PI/beta* p*nu/(denom.real()*denom.real()+denom.imag()*denom.imag());
    terms[k].second=nu;
  } 
  return terms;
}

void BCF_Decomposition_MeierTannor::print_info(std::ostream &ofs)const{
  ofs<<"DrudeLorentz: p="<<p<<" Omega="<<Omega<<" Gamma="<<Gamma<<" beta="<<beta<<std::endl;
}


void BCF_Decomposition_MeierTannor::setup(Parameters &param, const std::string &prefix){

  auto add_prefix = [&prefix] (const std::string &str){
    if(prefix=="")return str;
    return std::string(prefix+"_"+str);
  };

  //default units: ps^{-1}
  p=param.get_as_double(add_prefix("p"),1.);
  Omega=param.get_as_double(add_prefix("Omega"),1.);
  Gamma=param.get_as_double(add_prefix("Gamma"),1.);
  double temperature=readTemperature(param, prefix);
  beta=param.get_as_double(add_prefix("beta"),hbar_in_meV_ps/(kB_in_meV_by_K*temperature));
}


} //namespace
