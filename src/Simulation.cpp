#include "Simulation.hpp"
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
  void Simulation_Template<STATE, INFLUENCE>::
    run_return(Propagator &prop, const INFLUENCE &IF, 
           double ta, double dt, double te, Eigen::MatrixXcd &rho, 
           bool silent, bool use_symmetric_Trotter){

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
      ADM.propagate(prop, IF, t_in, dt, step, compressor, use_symmetric_Trotter);
      if(!silent)ADM.print_status(std::cout);
   
      Hamil=prop.get_Htot(t_in);
      set_results(step, t_in+dt, ADM.rho, &Hamil);

      if(!silent)std::cout<<std::endl;
    }
    rho=ADM.rho;
  }

template <typename STATE, typename INFLUENCE> 
  void Simulation_Template<STATE, INFLUENCE>::
  run_nobath(Propagator &prop, 
           double ta, double dt, double te, Eigen::MatrixXcd rho, 
           bool silent){

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


template class Simulation_Template<AugmentedDensityMatrix, InfluenceFunctional>;

template class Simulation_Template<AugmentedDensityMatrix_MPS, InfluenceFunctional_Vector>; 

}//namespace
