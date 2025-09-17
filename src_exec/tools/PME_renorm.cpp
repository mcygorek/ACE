#include "PolaronME.h"
//#include "Simulation.hpp"
#include "ReadTemperature.hpp"
#include "ACE.hpp"
#include "TimeGrid.hpp"

using namespace ACE;

int main(int args, char** argv){
  Parameters param(args, argv, true);

  std::string prefix=param.get_as_string("prefix", "Boson");

  double temperature=readTemperature(param, prefix);
  EnergyRange erange(param,prefix);
  double Emin=erange.E_min();
  double Emax=erange.E_max();
  if(Emax<=0.){
    std::cerr<<prefix<<"_E_max <= 0 !"<<std::endl;
    exit(1);
  }
  int Ndiscr=param.get_as_int("Ndiscr", 1e5);
  std::cout<<"N_discr: "<<Ndiscr<<std::endl;
  std::cout<<"E_max: "<<Emax<<std::endl;

  if(!SpectralDensity_Selector::is_specified(param, std::string(prefix+"_J"))){
    param.add_if_not_specified(std::string(prefix+"_J_type"), "QDPhonon");
  }
  RealFunctionPtr SD=SpectralDensity_Selector(param, std::string(prefix+"_J"));
 
  double renorm=PolaronME::renorm(SD, temperature, Emax, Ndiscr);
  std::cout<<"<B>: "<<renorm<<std::endl;
  double shift_offset=param.get_as_double("shift_offset",0);
  std::cout<<"shift_offset: "<<shift_offset<<std::endl;
  double PS=PolaronME::polaron_shift(SD, Emin, Emax, Ndiscr, shift_offset);
  std::cout<<"polaron shift: "<<PS<<" meV ("<<PS/hbar_in_meV_ps<<" ps^{-1})"<<std::endl;
 

  double Omega=param.get_as_double("Omega",0.);
  double delta=param.get_as_double("delta",0.);
  double eta=sqrt(delta*delta/hbar_in_meV_ps/hbar_in_meV_ps+Omega*Omega);
  double gammaPD_prefac=1.;
  if(Omega>1e-8)gammaPD_prefac=Omega*Omega/(eta*eta);

  double gammaPD=gammaPD_prefac*M_PI/2.*SD->f(eta)*coth(
                 beta_from_T(temperature)*hbar_in_meV_ps*eta/2.);
  std::cout<<"gammaPD: "<<gammaPD<<std::endl;

  double eta2=eta*renorm;
  double gammaPD_pol=M_PI/2.*SD->f(eta2)*coth(
                 beta_from_T(temperature)*hbar_in_meV_ps*eta2/2.);
  std::cout<<"gammaPD(polaron): "<<gammaPD_pol<<std::endl;

  
  std::string print_coherences=param.get_as_string("print_coherences");
  if(print_coherences!=""){
    TimeGrid tgrid(param);
    std::ofstream ofs(print_coherences.c_str());
    for(int n=0; n<=tgrid.n_tot; n++){
      double t=n*tgrid.dt;
//      std::complex<double> C=renorm*renorm*exp(PolaronME::phi(t, SD, temperature, Emax, Ndiscr));
      std::complex<double> C=exp(PolaronME::phi(t, SD, temperature, Emax, Ndiscr)-PolaronME::phi(0, SD, temperature, Emax, Ndiscr));
      ofs<<t<<" "<<C.real()<<" "<<C.imag()<<std::endl;
    }
  }

  return 0;
}
