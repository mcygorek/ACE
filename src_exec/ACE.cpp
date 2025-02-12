#include "PCH.hpp"
#include "GenericSimulation.hpp"
#include "ProcessTensorForwardList.hpp"
#include "TransferTensor.hpp"
#include "DummyException.hpp"
#include "ParametersScan.hpp"
#include "Timings.hpp"

using namespace ACE;

int main(int args, char** argv){
 try{
  Parameters param(args, argv, true);

  TimeGrid tgrid(param);

  int nr_scans=param.get_as_int("nr_scans",1);
  if(nr_scans<1){
    std::cerr<<"nr_scans < 1!"<<std::endl;
    exit(1);

  }else if(nr_scans==1){

    bool dont_propagate=param.get_as_bool("dont_propagate",false);
    if(dont_propagate){ param.override_param("outfile","/dev/null"); }
   
      GenericSimulation sim(param);

      std::cout<<"Setting up Process Tensor..."<<std::endl;
      Parameters paramPT=TransferTensor::get_paramPT(param);

      ProcessTensorForwardList PT;
      try{
        PT.setup(paramPT, sim.sysdim);
      }catch(DummyException &e){
        std::cerr<<"while creating ProcessTensorForwardList"<<std::endl;
        throw e;
      }
      PT.print_info();

      if(dont_propagate){
        std::cout<<"Not propagating ('dont_propagate' was set)"<<std::endl;
      }else{
        std::cout<<"Propagating..."<<std::endl;
        time_point time1=now();
        sim.run(tgrid, PT);
        time_point time2=now();
        std::cout<<"runtime for propagation: "<<time_diff(time2-time1)<<"ms"<<std::endl;
      }
  }else{ //nr_scans>1:

  
  //Check if dimensions are reasonable before calculating PT:
    {
      Parameters param2=ParametersScan(param,0);
      int consistent_dim=Output_Ops(param2).get_dim();
      Eigen::MatrixXcd initial_rho=InitialState(param2);
      FreePropagator fprop(param2);
      if(initial_rho.rows()!=consistent_dim){
         std::cerr<<"Mismatch in dimensions between initial system state and system propagator!\nPlease specify 'initial' and check whether the dimensions agree with 'add_Hamiltonian', 'add_Lindblad', ..."<<std::endl;
        exit(1);
      }
      fprop.set_dim(consistent_dim, "Mismatch in dimensions between observables system propagator!\nPlease specify 'add_Output' and check whether the dimensions agree with 'add_Hamiltonian', 'add_Lindblad', ...");


      std::cout<<"Setting up process tensor..."<<std::endl;
      ProcessTensorForwardList PT(param, consistent_dim);
      PT.print_info();

      for(int scan=0; scan<nr_scans; scan++){
        std::cout<<std::endl;
        std::cout<<"###################################"<<std::endl;
        std::cout<<"Scan "<<scan<<"/"<<nr_scans<<std::endl;
        std::cout<<"###################################"<<std::endl;
        Parameters param2=ParametersScan(param,scan);
        GenericSimulation sim(param2);
        sim.run(tgrid, PT);
      } 
    }
  }
 }catch (DummyException &e){
  return 1;
 }
#ifdef EIGEN_USE_MKL_ALL
  mkl_free_buffers();
#endif

  return 0;
}

