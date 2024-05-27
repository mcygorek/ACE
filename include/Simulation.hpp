#ifndef SIMULATION_DEFINED_H
#define SIMULATION_DEFINED_H

#include "Parameters.hpp"
#include "Operators.hpp"
#include "FreePropagator.hpp"
#include "InfluenceFunctional.hpp"
#include "ADM.hpp"
#include <fstream>  
#include "Simulation_Results.hpp"
#include "FT_Output.hpp"
#include "InitialState.hpp"
#include "SpectralDensity_Selector.hpp"
#include "RankCompressor_Selector.hpp"

#include "InfluenceFunctional_Vector.hpp"
#include "ADM_MPS.hpp"

namespace ACE{

template <typename STATE, typename INFLUENCE>
class Simulation_Template{
public:

  Simulation_Results results;
  ///defines which operator averages are to be printed
  Output_Ops output_Op;
  FT_Output FTO;

  Smart_Ptr<RankCompressor> compressor;

  inline void add_output_Op(const Eigen::MatrixXcd &op){
    output_Op.add(op);
  }
  inline void add_output_proj(int i){
    output_Op.add_proj(i);
  }
  inline void setup_output(Parameters &param){
    output_Op.setup(param);
    FTO.setup(param);
  }

  inline void set_results(int step, double t, const Eigen::MatrixXcd &rho, const Eigen::MatrixXcd *H=NULL){
    results.set(step, t, output_Op, rho, H);
  }
  inline void set_compressor(Smart_Ptr<RankCompressor>  &compr_ptr){
    compressor=compr_ptr;
  }
  inline void set_compressor(RankCompressor *compr_ptr){
    compressor=compr_ptr;
  }

  void run_return(Propagator &prop, const INFLUENCE &IF, 
           double ta, double dt, double te, Eigen::MatrixXcd &rho, 
           bool silent, bool use_symmetric_Trotter);
  
  inline void run(Propagator &prop, const INFLUENCE &IF, 
           double ta, double dt, double te, Eigen::MatrixXcd rho, 
           bool silent, bool use_symmetric_Trotter){
     run_return(prop, IF, ta, dt, te, rho, silent, use_symmetric_Trotter);
  }
  inline void run_silent(Propagator &prop, const INFLUENCE &IF, 
                  double ta, double dt, double te, const Eigen::MatrixXcd &rho,
                  bool use_symmetric_Trotter){
    run(prop, IF, ta, dt, te, rho, true, use_symmetric_Trotter);
  }
  void run_nobath(Propagator &prop, 
           double ta, double dt, double te, Eigen::MatrixXcd rho, 
           bool silent=false);
  
  inline void print_results(const std::string &fname)const{
    results.print(fname);
    FTO.print(results);
  }

};


typedef Simulation_Template<AugmentedDensityMatrix, InfluenceFunctional> Simulation;

typedef Simulation_Template<AugmentedDensityMatrix_MPS, InfluenceFunctional_Vector> Simulation_TEMPO; 

}//namespace
#endif
