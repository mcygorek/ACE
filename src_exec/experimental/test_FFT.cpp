#include "Parameters.hpp"
#include "discreteFT.hpp"
#include "Simulation_Results.hpp"
#include <iostream>
#include <cstdlib>

using namespace ACE;

int main(int args, char** argv){
  Parameters param(args, argv, true);

  int N=param.get_as_int("N", 1024);
  double t_min=param.get_as_double("t_min",  0.);
  double t_max=param.get_as_double("t_max", 10.);
  double omega_min=param.get_as_double("omega_min",  0.);
  double omega_max=param.get_as_double("omega_max", 10.);
  double domega_min=param.get_as_double("domega_min",  0.);
  int sign=param.get_as_int("sign", -1);

  double center=param.get_as_double("center", 0); //t_max/2.);
  double sigma=param.get_as_double("sigma", 1.);
  std::string print_input=param.get_as_string("print_input","");
  std::string print_output=param.get_as_string_check("print_output");
  int integrate_mode=param.get_as_int("integrate_mode",0);


  std::cerr<<"N="<<N<<std::endl;
  std::cerr<<"t_max="<<t_max<<std::endl;
  double dt=t_max/N;

/*
  int lN=(log(N)/log(2)+0.5);
  if(N!=pow(2,lN)){ 
    std::cerr<<"N="<<N<<" is not a power of 2!"<<std::endl;
    exit(1);
  }
  std::cout<<"N=2^"<<lN<<std::endl;
*/
  Simulation_Results in(N);  
  for(size_t i=0; i<in.size(); i++){
    double t=t_min+i*dt;
    double x=(t-center)/sigma;

    std::complex<double> c=exp(-0.5*x*x);
//    std::complex<double> c=cos(10.*t);

    in[i].first=t;
    in[i].second=std::vector<std::complex<double> > (1, c);
  }

  if(print_input!=""){
    in.print(print_input);
  }

  Simulation_Results out=resultsFFT(in, 0, sign, integrate_mode);

  out.print(print_output);

  return 0;
}
