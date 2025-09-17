#ifndef ACE_OUTPUT_PRINTER_DEFINED_H
#define ACE_OUTPUT_PRINTER_DEFINED_H

#include "Output_Ops.hpp"
#include "Which_Env_Ops.hpp"
#include "Parameters.hpp"  
#include "ProcessTensorElement.hpp"
#include <fstream>
#include <memory>

namespace ACE{

class OutputPrinter{
public:  

  bool print_timestep;
  std::unique_ptr<std::ofstream> ofs;

  bool full_densmat; //extracts full density matrix; overrides Output_Ops
  Output_Ops output_Op;
  Which_Env_Ops_List which_env_ops;

//Print occupation of instantaneous eigenstates:
  std::unique_ptr<std::ofstream> ofs_eigenstates;
  std::vector<size_t> eigenstate_components;

//Storage for reduced density matrices at last few time steps:
  std::vector<Eigen::VectorXcd> rho_t;
  std::vector<double> rho_times;
  int start_extract; 
  bool do_extract;


  virtual void clear();
  virtual void setup(Parameters & param, int setdim=-1);
  virtual void set_stream(const std::string &outfile, int precision=-1);

  virtual void print(int n, double t, const Eigen::VectorXcd & rho_reduced,
                 const std::vector<std::complex<double> > & env_reduced = 
                                      std::vector<std::complex<double> >() );

  virtual void print_eigenstate_occupations(double t, const Eigen::MatrixXcd &H, const Eigen::MatrixXcd &prop);

  virtual void finish();

  virtual std::pair<Eigen::VectorXd,Eigen::MatrixXcd> extract()const;
//  inline const std::vector<Eigen::VectorXcd> & extract()const{
//    return rho_t;
//  }

  OutputPrinter(const std::string &outfile, const std::vector<Eigen::MatrixXcd> & list){
    clear();
    set_stream(outfile);
    output_Op=Output_Ops(list);
  }
  OutputPrinter(Parameters & param, int setdim=-1){
    setup(param, setdim);
  }
  OutputPrinter(){
    clear();
  }
};
}//namespace
#endif
