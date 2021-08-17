#ifndef SIMULATION_DEFINED_H
#define SIMULATION_DEFINED_H

#include "Parameters.h"
#include "Operators.h"
#include "FreePropagator.h"
#include "InfluenceFunctional.h"
#include "ADM.h"
#include <fstream>  
#include "Simulation_Results.h"
#include "FT_Output.h"
#include "InitialState.h"
#include "SpectralDensity_Selector.h"
#include "RankCompressor_Selector.h"

#include "InfluenceFunctional_Vector.h"
#include "ADM_MPS.h"

template <typename STATE, typename INFLUENCE>
class Simulation_Template{
public:

  Simulation_Results results;
  ///defines which operator averages are to be printed
  Output_Ops output_Op;
  FT_Output FTO;

  Smart_Ptr<RankCompressor> compressor;

  void add_output_Op(const Eigen::MatrixXcd &op){
    output_Op.add(op);
  }
  void add_output_proj(int i){
    output_Op.add_proj(i);
  }
  void setup_output(Parameters &param){
    output_Op.setup(param);
    FTO.setup(param);
  }

  void set_results(int step, double t, const Eigen::MatrixXcd &rho, const Eigen::MatrixXcd *H=NULL){
    results.set(step, t, output_Op, rho, H);
  }
  void set_compressor(Smart_Ptr<RankCompressor>  &compr_ptr){
    compressor=compr_ptr;
  }
  void set_compressor(RankCompressor *compr_ptr){
    compressor=compr_ptr;
  }

  void run_return(Propagator &prop, const INFLUENCE &IF, 
           double ta, double dt, double te, Eigen::MatrixXcd &rho, 
           bool silent=false){

    prop.check_dimensions();
    if(rho.cols()!=rho.rows()){
      std::cerr<<"Simulation::calculate: rho.cols()!=rho.rows()!"<<std::endl;
      exit(1);
    }
    if(rho.rows()!=prop.get_dim()){
      std::cerr<<"Simulation::calculate: rho.rows()!=prop.get_dim()!"<<std::endl;
      exit(1);
    }

//    int dim=prop.get_dim();
    size_t Nsteps=(te-ta)/dt;


    STATE ADM(IF.get_n_max(), IF.get_Ngrps(), rho);

    results.clear(); 
    results.resize(Nsteps+1);
    Eigen::MatrixXcd Hamil=prop.get_Htot(ta);
    set_results(0, ta, rho, &Hamil);
    for(size_t step=1; step<=Nsteps; step++){
      if(!silent)std::cout<<"step: "<<step<<"/"<<Nsteps<<" "<<std::flush;

      double t_in=ta+(step-1)*dt;
      ADM.propagate(prop, IF, t_in, dt, step, compressor);
      ADM.print_status(std::cout);
   
      Hamil=prop.get_Htot(t_in);
      set_results(step, t_in+dt, ADM.rho, &Hamil);

      if(!silent)std::cout<<std::endl;
    }
    rho=ADM.rho;
  }
  void run(Propagator &prop, const INFLUENCE &IF, 
           double ta, double dt, double te, Eigen::MatrixXcd rho, 
           bool silent=false){
     run_return(prop, IF, ta, dt, te, rho, silent);
  }
  void run_silent(Propagator &prop, const INFLUENCE &IF, 
                  double ta, double dt, double te, const Eigen::MatrixXcd &rho){
    run(prop, IF, ta, dt, te, rho, true);
  }
  void run_nobath(Propagator &prop, 
           double ta, double dt, double te, Eigen::MatrixXcd rho, 
           bool silent=false){

    prop.check_dimensions();
    if(rho.cols()!=rho.rows()){
      std::cerr<<"Simulation::calculate: rho.cols()!=rho.rows()!"<<std::endl;
      exit(1);
    }
    if(rho.rows()!=prop.get_dim()){
      std::cerr<<"Simulation::calculate: rho.rows()!=prop.get_dim()!"<<std::endl;
      exit(1);
    }

    int dim=prop.get_dim();
    size_t Nsteps=(te-ta)/dt;
std::cout<<"RUN_NOBATH!"<<std::endl;

    results.clear(); 
    results.resize(Nsteps+1);
    set_results(0, ta, rho);
    for(size_t step=1; step<=Nsteps; step++){
      if(!silent)std::cout<<"step: "<<step<<"/"<<Nsteps<<std::endl;

      double t_in=ta+(step-1)*dt;
      prop.update(t_in,dt);

      Eigen::MatrixXcd rho_tmp=Eigen::MatrixXcd::Zero(dim, dim);
      for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
          for(int k=0; k<dim; k++){
            for(int l=0; l<dim; l++){
              rho_tmp(i,j) +=  prop.M(i*dim+j,k*dim+l) * rho(k,l);
            }
          }
        }
      }
   
      rho=rho_tmp;
      set_results(step, t_in+dt, rho);
    }

  }

  void print_results(const std::string &fname)const{
    results.print(fname);
    FTO.print(results);
  }

};


typedef Simulation_Template<AugmentedDensityMatrix, InfluenceFunctional> Simulation;

typedef Simulation_Template<AugmentedDensityMatrix_MPS, InfluenceFunctional_Vector> Simulation_TEMPO; 


#endif
