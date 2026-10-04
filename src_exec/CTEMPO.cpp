#include "Simulation_CTEMPO.hpp"
#include "Timings.hpp"


using namespace ACE;

int main(int args, char** argv){
 try{
  Parameters param(args, argv, true);

  if(param.is_specified("N_sites")){
    std::cerr<<"CTEMPO called with parameter 'N_sites'. Did you want to run CCTEMPO?"<<std::endl;
    throw DummyException();
  }

  TimeGrid tgrid(param);
  InitialState initial(param);
  FreePropagator fprop(param);
  if(!fprop.dim_fixed){
    fprop.set_dim(initial.rho.rows());
  }else{
    if(fprop.get_dim()!=initial.rho.rows()){
      std::cerr<<"Mismatch in dimensions between initial state and propagator!"<<std::endl;
      exit(1);
    }
  }

  std::string prefix=param.get_as_string("prefix","Boson");
  DiagBB diagBB; diagBB.setup(param,prefix);  
  OutputPrinter printer(param);
  TruncationLayout trunc_layout(param);
  Simulation_CTEMPO sim(param);

  time_point time1=now();
  sim.run(fprop, diagBB, initial.rho, tgrid, printer, trunc_layout);
  time_point time2=now();
  std::cout<<"CTEMPO runtime: "<<time_diff(time2-time1)<<"ms"<<std::endl;

 }catch (DummyException &e){
  return 1;
 }

#ifdef EIGEN_USE_MKL_ALL
  mkl_free_buffers();
#endif

  return 0;
}

