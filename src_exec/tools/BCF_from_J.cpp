#include "Parameters.hpp"
#include "EnergyRange.hpp"
#include "Constants.hpp"
#include "ReadTemperature.hpp"
#include "SpectralDensity_Selector.hpp"

using namespace ACE;

int main(int args, char** argv){

  Parameters param(args, argv, true);
  double dt=param.get_as_double("dt",0.01);
  double te=param.get_as_double("te",20);
  std::string print_BCF=param.get_as_string("print_BCF", "BCF.dat");
  std::string Gaussian_prefix=param.get_as_string("Gaussian_prefix", "Boson");

  RealFunctionPtr SD=SpectralDensity_Selector(param, Gaussian_prefix+"_J");

  EnergyRange Erange(param, Gaussian_prefix);
  double wa=Erange.omega_min();
  double we=Erange.omega_max();
  int Ndiscr=param.get_as_size_t("Ndiscr", 1e5);
  if(Ndiscr%2==1){
    std::cerr<<"Please set 'Ndiscr' to an even number"<<std::endl;
    exit(1);
  }
  double dw=(we-wa)/Ndiscr;
  if(dw<1e-16){
    std::cerr<<"Please set 'Boson_omega_max' to a value larger than 'Boson_omega_min'!"<<std::endl;
    exit(1);
  }
  double temperature=readTemperature(param,Gaussian_prefix);
  double beta=hbar_in_meV_ps/(kB_in_meV_by_K*temperature);

  std::cout<<"wa="<<wa<<" we="<<we<<" dw="<<dw<<std::endl;
  std::cout<<"beta="<<beta<<std::endl;


  std::ofstream ofs(print_BCF.c_str());
  int NT=(te/dt+0.5);
  for(int l=0; l<=NT; l++){
    double t=l*dt;
    std::complex<double> c=0;

    for(int k=0; k<Ndiscr; k++){
      double w=wa+(k+0.5)*dw;
      std::complex<double> contr=dw*SD->f(w)*exp(std::complex<double>(0.,-w*t));
      if(temperature>1e-10){
        contr.real(contr.real()*coth(beta*w/2.));
      }
      c+=contr;
    }
    ofs<<t<<" "<<c.real()<<" "<<c.imag()<<std::endl;
  }

  std::cout<<"Bath correlation function written to file '"<<print_BCF<<"'."<<std::endl;
}
