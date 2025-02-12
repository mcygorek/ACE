#include "BCF_Decomposition_DrudeLorentz.hpp"
#include "Parameters.hpp"

using namespace ACE;

int main(int args, char** argv){

  Parameters param(args, argv, true);
  int Nmats=param.get_as_size_t("N_Matsubara",10);
  double dt=param.get_as_double("dt",0.01);
  double te=param.get_as_double("te",20);
  std::string print_BCF=param.get_as_string("print_BCF", "BCF.dat");
  std::string print_J=param.get_as_string("print_J", "J.dat");
  double print_J_wa=param.get_as_double("print_J", 0, 0, 1);
  double print_J_we=param.get_as_double("print_J", 100, 0, 2);
  int print_J_Ndiscr=param.get_as_size_t("print_J", 10000, 0, 3);

  BCF_Decomposition_DrudeLorentz bcfd(param);

  bcfd.print_info();
  bcfd.print_all(Nmats);
  bcfd.print_BCF(print_BCF, dt, te, Nmats);
  bcfd.print_J(print_J, print_J_wa, print_J_we, print_J_Ndiscr);
}
