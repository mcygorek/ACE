#ifndef ACE_OUTPUT_PRINTER_CCTEMPO_DEFINED_H
#define ACE_OUTPUT_PRINTER_CCTEMPO_DEFINED_H

#include "Parameters.hpp"
#include "FreePropagator.hpp"
#include "MPO_tools.hpp"
#include "Comb.hpp"
#include <iomanip>

namespace ACE{
struct OutputPrinter_CCTEMPO{
  std::unique_ptr<std::ofstream> ofs;
  //print order: First OneBody in order of sites, then collective, then coeff_list
  std::vector<std::pair<int, Eigen::VectorXcd> > OneBody;
  
  std::vector<MPO> collective;

  //list of explicit coefficients to be printed
  std::vector<std::vector<int> > coeff_list;

  //Print occupation of instantaneous eigenstates:
  std::unique_ptr<std::ofstream> ofs_eigenstates;
  std::vector<size_t> eigenstate_components;



  void set_stream(const std::string &outfile, int precision);
  void clear();

  void print(int n, double t, const MPO & state, const std::vector<Eigen::VectorXcd> &closure);
  void print(int n, double t, Comb &comb);
  
  void print_eigenstate_occupations(double t, const Eigen::MatrixXcd & H, const MPO & state, const std::vector<Eigen::VectorXcd> & closure);
  void print_eigenstate_occupations(double t, const Eigen::MatrixXcd & H, Comb &comb);



  void setup(Parameters & param);
  OutputPrinter_CCTEMPO(Parameters & param){
    setup(param);
  }
  OutputPrinter_CCTEMPO(){
    clear();
  }
};

}//namespace
#endif
