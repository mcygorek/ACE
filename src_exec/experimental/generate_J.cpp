#include "SpectralDensity_Selector.hpp"
#include "MPG_Discretization.hpp"
#include "Parameters.hpp"
#include "Reader.hpp"

using namespace ACE;

int main(int args, char** argv){
  Parameters param(args, argv, true);
  
  std::string prefix=param.get_as_string("prefix","");

  SpectralDensity_Selector SD_select(param,prefix);
  RealFunctionPtr SD=SD_select;


  EnergyRange range(param);
  if(fabs(range.E_range())<1e-12){
    std::cerr<<"Discretization domain too small! Please specify 'E_max' and/or 'E_min'!"<<std::endl; 
    exit(1);
  }
  
  int N=param.get_as_int(add_prefix(prefix,"N"), 1e5);

  std::string outfile=param.get_as_string("outfile","J.dat");
  SD->print(outfile, range.omega_min(), range.omega_max(), N);
 
  std::cout<<"Spectral density written to file '"<<outfile<<"'."<<std::endl;
  
  return 0;
}
