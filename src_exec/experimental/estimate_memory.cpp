#include "Parameters.hpp"
#include <iostream>
#include <cstdlib>
#include "DiagBB.hpp"
#include "TimeGrid.hpp"
#include "ProcessTensorBuffer.hpp"

using namespace ACE;

int main(int args, char** argv){
  Parameters param(args, argv, true);

  std::string prefix=param.get_as_string("prefix","Boson");
  param.complain_if_not_specified("dt");
  param.complain_if_not_specified("te");
  param.complain_if_not_specified("threshold");

  TimeGrid tgrid(param);
  TruncatedSVD trunc(param);
  trunc.print_info(); std::cout<<std::endl;
  DiagBB diagBB(param, prefix);

  ProcessTensorBuffer PTB;
  PTB.set_from_DiagBB_single_line_auto(diagBB, tgrid.dt, tgrid.n_tot, -2, trunc, true, false);
  int n_mem=PTB.get_n_tot();

  std::cout<<"Estimated memory time: t_mem="<<n_mem*tgrid.dt<<" (n_mem="<<n_mem<<") of total time "<<tgrid.get_t_tot()<<" (n_tot="<<tgrid.n_tot<<")"<<std::endl;

  return 0;
}
