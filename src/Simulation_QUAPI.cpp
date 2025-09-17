#include "Simulation_QUAPI.hpp"
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


    STATE ADM(IF.get_n_max()-1, IF.get_Ngrps(), rho);

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
    run(Propagator &prop, const INFLUENCE &IF, const Eigen::MatrixXcd &rho, 
        const TimeGrid &tgrid, OutputPrinter &printer, 
        bool use_symmetric_Trotter, bool silent){

    prop.check_dimensions();
    if(rho.cols()!=rho.rows()){
      std::cerr<<"Simulation::calculate: rho.cols()!=rho.rows()!"<<std::endl;
      exit(1);
    }
    if(rho.rows()!=prop.get_dim()){
      std::cerr<<"Simulation::calculate: rho.rows()!=prop.get_dim()!"<<std::endl;
      exit(1);
    }

    size_t Nsteps=(tgrid.get_t_tot()-tgrid.ta)/tgrid.dt;

    STATE ADM(IF.get_n_max()-1, IF.get_Ngrps(), rho);

/*
if(printer.ofs->is_open()){
std::cout<<"printer.ofs->is_open()=true"<<std::endl;
}else{
std::cout<<"printer.ofs->is_open()=false"<<std::endl;
}
*printer.ofs<<"#test line"<<std::endl;
std::cout<<"test line printed to outfile"<<std::endl;
if(printer.do_extract){
std::cout<<"printer.do_extract=true"<<std::endl;
}else{
std::cout<<"printer.do_extract=false"<<std::endl;
}
std::cout<<"printer.output_Op.size()="<<printer.output_Op.size()<<std::endl;
*/


    printer.print(0, tgrid.get_t(0), H_Matrix_to_L_Vector(rho));
    printer.print_eigenstate_occupations(tgrid.get_t(0), prop.get_Htot(tgrid.get_t(0)), rho);

    for(size_t step=1; step<=Nsteps; step++){
      if(!silent)std::cout<<"step: "<<step<<"/"<<Nsteps<<" "<<std::flush;
      double t_in=tgrid.ta+(step-1)*tgrid.dt;
      ADM.propagate(prop, IF, t_in, tgrid.dt, step, compressor, use_symmetric_Trotter);
      if(!silent)ADM.print_status(std::cout);
   
      printer.print(step, t_in+tgrid.dt, H_Matrix_to_L_Vector(ADM.rho));
      printer.print_eigenstate_occupations(t_in+tgrid.dt, prop.get_Htot(t_in+tgrid.dt), ADM.rho);


      if(!silent)std::cout<<std::endl;
    }
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
