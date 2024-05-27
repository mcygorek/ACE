#include "BCF_Decomposition_DrudeLorentz.hpp"
#include "ReadTemperature.hpp"
#include "Constants.hpp"
#include <fstream>


namespace ACE{

double BCF_Decomposition_DrudeLorentz::J(double omega)const{
  return (2./M_PI) * lambda*gamma*omega/(gamma*gamma+omega*omega);
}


std::vector<BCF_Decomposition::Term> BCF_Decomposition_DrudeLorentz::get_nonMatsubara()const{
  std::vector<BCF_Decomposition::Term> terms(1);
  
  terms[0].first=lambda*gamma*std::complex<double>(1./tan(beta*gamma/2.),-1.);
  terms[0].second=gamma;
  
  return terms;
}
std::vector<BCF_Decomposition::Term> BCF_Decomposition_DrudeLorentz::get_Matsubara(int N)const{
  std::vector<BCF_Decomposition::Term> terms(N);
  for(int k=0; k<N; k++){
    double nu=2.*M_PI*(k+1)/beta;
    terms[k].first=4.*lambda*gamma*nu/((nu*nu-gamma*gamma)*beta);
    terms[k].second=nu;
  } 
  return terms;
}

bool BCF_Decomposition_DrudeLorentz::use_terminator()const{
  return true;
}
std::complex<double> BCF_Decomposition_DrudeLorentz::get_terminator(int Nk)const{
  std::complex<double> ret(2.*lambda/(beta*gamma),-lambda);

  std::vector<Term> terms=get_nonMatsubara();
  for(const auto &t :terms){
    ret-=t.first/t.second;
  }
  int nr_nM=terms.size();

  terms=get_Matsubara(Nk-nr_nM);
  for(const auto &t :terms){
    ret-=t.first/t.second;
  }
  
  return 0.5*ret;
}

void BCF_Decomposition_DrudeLorentz::print_info(std::ostream &ofs)const{
  ofs<<"DrudeLorentz: lambda="<<lambda<<" gamma="<<gamma<<" beta="<<beta<<std::endl;
}


void BCF_Decomposition_DrudeLorentz::setup(Parameters &param, const std::string &prefix){

  auto add_prefix = [&prefix] (const std::string &str){
    if(prefix=="")return str;
    return std::string(prefix+"_"+str);
  };

  //default units: ps^{-1}
  lambda=param.get_as_double(add_prefix("lambda"),1.);
  gamma=param.get_as_double(add_prefix("gamma"),1.);
  double temperature=readTemperature(param, prefix);
  beta=param.get_as_double(add_prefix("beta"),hbar_in_meV_ps/(kB_in_meV_by_K*temperature));
}


} //namespace
