#include "Parameters.hpp"
#include "ReaderBasics.hpp"
#include "TimeGrid.hpp"
#include <fstream>

using namespace ACE;

int main(int args, char **argv){
  Parameters param(args, argv);

  std::string outfile=param.get_as_string_check("outfile");
  double g=param.get_as_double_check("g");
  double gamma=param.get_as_double_check("gamma");
  double omega=param.get_as_double_check("omega");
  TimeGrid tgrid(param);
  double dt=tgrid.dt;

  std::ofstream ofs(outfile.c_str());
  std::complex<double> k(gamma, omega);
  if(abs(k)<1e-6){
    std::cerr<<"|gamma+i*omega|<1e-6!"<<std::endl;
    exit(1); 
  }
  std::complex<double> c=g*g*(dt/k+(exp(-k*dt)-1.)/(k*k));
  ofs<<0<<" "<<c.real()<<" "<<c.imag()<<std::endl;
  for(int n=1; n<tgrid.n_tot; n++){
    double t=n*dt;
    c=g*g*(exp(k*dt)+exp(-k*dt)-2.)/(k*k)*exp(-k*t);
    ofs<<t<<" "<<c.real()<<" "<<c.imag()<<std::endl;
  }

  return 0;
}
