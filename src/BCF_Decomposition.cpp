#include "BCF_Decomposition.hpp"
#include <fstream>

namespace ACE{

double BCF_Decomposition::J(double omega)const{
  return 0;
}

std::vector<BCF_Decomposition::Term> BCF_Decomposition::get_all(int N_matsubara)const{
  std::vector<BCF_Decomposition::Term> terms=get_nonMatsubara();
  std::vector<BCF_Decomposition::Term> terms2=get_Matsubara(N_matsubara);
  terms.insert(terms.end(), terms2.begin(), terms2.end());
  return terms;
}

void BCF_Decomposition::print_all(int N_matsubara, std::ostream &ofs)const{
  std::vector<BCF_Decomposition::Term> terms=get_all(N_matsubara);
  for(const BCF_Decomposition::Term &t: terms){
    ofs<<t.first<<" "<<t.second<<std::endl;
  }
}
void BCF_Decomposition::print_BCF(const std::string &fname, double dt, double te, int N_matsubara)const{
  int NT=(te/dt+0.5);
  std::ofstream ofs(fname.c_str());

  std::vector<BCF_Decomposition::Term> terms=get_all(N_matsubara);
  for(int l=0; l<=NT; l++){
    double t=l*dt;
    std::complex<double> bfc=0.;
    for(const BCF_Decomposition::Term &term: terms){
      bfc+=term.first*exp(-term.second*t);
    }
    ofs<<t<<" "<<bfc.real()<<" "<<bfc.imag()<<std::endl;
  }
} 

void BCF_Decomposition::print_J(const std::string &str, double xa, double xb, int Ndiscr)const{
  std::ofstream ofs(str.c_str());
  double h=(xb-xa)/Ndiscr;
  for(int i=0; i<=Ndiscr; i++){
    double x=xa+i*h;
    ofs<<x<<" "<<J(x)<<std::endl;
  }
}


bool BCF_Decomposition::use_terminator()const{
  return false;
}
std::complex<double> BCF_Decomposition::get_terminator(int Nk)const{
 return 0.;
}

}
