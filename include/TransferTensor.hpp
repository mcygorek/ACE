#ifndef ACE_TRANSFER_TENSOR_DEFINED_H
#define ACE_TRANSFER_TENSOR_DEFINED_H

#include "OutputPrinter.hpp"
#include "Simulation_PT.hpp"
#include "ProcessTensorForwardList.hpp"
#include "TimeGrid.hpp"
#include <iostream>
#include <memory>

namespace ACE {

class TransferTensor{
public:

  std::vector<Eigen::MatrixXcd> T;
  Eigen::MatrixXcd LT; bool use_LT; //single long-time propagator

  static std::vector<Eigen::MatrixXcd> calculate_E(Propagator &prop, ProcessTensorForwardList &PT, Simulation_PT &sim, const TimeGrid &tgrid, int verbosity=0);

  void set_T_from_E(const std::vector<Eigen::MatrixXcd> &E, int n_max=-1);

  void set_LT_from_E(const std::vector<Eigen::MatrixXcd> &E, int verbosity=0);


  void calculate(Propagator &prop, ProcessTensorForwardList &PT, Simulation_PT &sim, const TimeGrid &tgrid);


  void propagate_initial(std::vector<Eigen::VectorXcd> &rho_t, const TimeGrid &tgrid, OutputPrinter &printer)const;
  void propagate_LT(std::vector<Eigen::VectorXcd> &rho_t, const TimeGrid &tgrid, OutputPrinter &printer)const;
  void propagate_bulk(std::vector<Eigen::VectorXcd> &rho_t, const TimeGrid &tgrid, OutputPrinter &printer)const;
  void propagate(const Eigen::MatrixXcd &initial_rho, const TimeGrid &tgrid, OutputPrinter &printer)const;
 
 
  //read and write not implemented yet:
  void read(const std::string &fname); 
  void write(const std::string &fname)const; 

  void print_norms(std::ostream &os=std::cout, double dt=1.)const;
  void print_norms(const std::string &fname, double dt=1.)const;

  static bool get_use_TT(Parameters &param);
  static int get_TT_n_from(Parameters &param);
  static int get_TT_n_mem(Parameters &param, bool do_complain);
//if TTs are use, PT only have to be calculated up to te'<te; produce corresponding paramPT for PT simulation:
  static Parameters get_paramPT(Parameters &param);


  TransferTensor(): use_LT(false){}
  TransferTensor(const std::string &fname): use_LT(false){
    read(fname);
  }
  TransferTensor(Propagator &prop, ProcessTensorForwardList &PT, Simulation_PT &sim, const TimeGrid &tgrid): use_LT(false){
    calculate(prop, PT, sim, tgrid);
  }

};
}//namespace


#endif
