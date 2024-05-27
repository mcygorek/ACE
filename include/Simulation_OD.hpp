#ifndef SIMULATION_OD_DEFINED_H
#define SIMULATION_OD_DEFINED_H

#include <memory>
#include "IF_from_Parameters.hpp"
#include "Propagator.hpp"
#include "Simulation_Results.hpp"
//#include "Simulation.h"
#include "InitialState.hpp"
#include "FT_Output.hpp"
#include "Which_Env_Ops.hpp"
//#include "Cavityfy.h"

#include "ModePropagatorGenerator.hpp"
//#include "Pulse_Printer.h"

namespace ACE{

class Simulation_OD{
public:

  Eigen::MatrixXcd rho;
  bool print_timestep;
  bool use_symmetric_Trotter;

  Simulation_Results results;
  Output_Ops output_Op;
  Which_Env_Ops_List which_env_ops;
  FT_Output FTO;
  

  void add_output_Op(const Eigen::MatrixXcd &op);
  
  void setup_output(Parameters &param);

  void run_nobath_dt0_nmax(Propagator &prop,
                           double ta, double dt, double dt0, int n_max, 
                           const Eigen::MatrixXcd &initial_rho);

  void run_nobath(Propagator &prop,
    double ta, double dt, double te, const Eigen::MatrixXcd &initial_rho);
    
  void run_dt0_nmax(Propagator &prop, 
                    std::vector<std::shared_ptr<IF_OD_Abstract> > & IFV,
                    double ta, double dt, double dt0, int n_max, 
                    const Eigen::MatrixXcd &initial_rho);
 
  void run_dt0_nmax_sym_multiPT(Propagator &prop,
                    std::vector<std::shared_ptr<IF_OD_Abstract> > & IFV,
                    double ta, double dt, double dt0, int n_max,
                    const Eigen::MatrixXcd &initial_rho);

  void run_sym_multiPT(Propagator &prop, 
    std::vector<std::shared_ptr<IF_OD_Abstract> > & IFV,
    const TimeGrid &tgrid, const Eigen::MatrixXcd &initial_rho);

  void run_sym_multiPT(Propagator &prop,
    std::vector<std::shared_ptr<InfluenceFunctional_OD> > & IFV,
    const TimeGrid &tgrid, const Eigen::MatrixXcd &initial_rho);


  void run(Propagator &prop, 
    std::vector<std::shared_ptr<IF_OD_Abstract> > & IFV,
    double ta, double dt, double te, const Eigen::MatrixXcd &initial_rho);
    
  void run(Propagator &prop, 
    std::vector<std::shared_ptr<InfluenceFunctional_OD> > & IFV,
    double ta, double dt, double te, const Eigen::MatrixXcd &initial_rho);
    
  void run(Propagator &prop, 
    std::vector<std::shared_ptr<IF_OD_Abstract> > & IFV,
    const TimeGrid &tgrid, const Eigen::MatrixXcd &initial_rho);
    
  void run(Propagator &prop,
    std::vector<std::shared_ptr<InfluenceFunctional_OD> > & IFV,
    const TimeGrid &tgrid, const Eigen::MatrixXcd &initial_rho);

  void run(Propagator &prop, std::shared_ptr<IF_OD_Abstract> IF,
    double ta, double dt, double te, const Eigen::MatrixXcd &initial_rho);
  
  void print_results(const std::string &fname)const;

  void initialize();
  
  inline Simulation_OD(Propagator &prop, std::shared_ptr<IF_OD_Abstract> IF, 
    double ta, double dt, double te, const Eigen::MatrixXcd &initial_rho){
     
    initialize();
    run(prop, IF, ta, dt, te, initial_rho);
  }
  inline Simulation_OD(){
    initialize();
  }
};

}//namespace
#endif
