#include "ProcessTensorRepeat.hpp"
#include "TempFileName.hpp"
#include "Parameters.hpp"
#include "GenericSimulation.hpp"
#include "ProcessTensorForwardList.hpp"
#include "DummyException.hpp"

using namespace ACE;
int main(int args, char** argv){
 try{
  Parameters param(args, argv, true);

  std::string read_PT = param.get_as_string_check("read_PT");
 
  double dict_zero = param.get_as_double("dict_zero",-1.);
  bool flip_last = param.get_as_bool("flip_last",true);
  bool sweep_initial = param.get_as_bool("sweep_initial",true);
  bool start_with_forward = param.get_as_bool("start_with_forward",false);
  int n_forward_sweep = param.get_as_size_t("n_forward_sweep",0);
  int n_sweep = param.get_as_size_t("n_sweep",(n_forward_sweep?1:0));

  std::cout<<"File: "<<read_PT<<std::endl;
  std::cout<<"n_sweep="<<n_sweep;
  if(flip_last){std::cout<<" flip_last=true";}else{std::cout<<" flip_last=false";}
  if(sweep_initial){std::cout<<" sweep_initial=true";}else{std::cout<<" sweep_initial=false";}
  std::cout<<std::endl;
  
  if(n_sweep>1 && !flip_last){
    std::cerr<<"Don't set flip_last=true when also using n_sweeps>1"<<std::endl<<std::endl;
    throw DummyException();
  }
  
//  int sweep_n = param.get_as_int("sweep_n",param.get_as_int("intermediate_sweep_n",0));
  int verbosity = param.get_as_int("verbosity", 2); 

  TruncationLayout trunc(param);
  trunc.print_info(); std::cout<<std::endl;
  TruncatedSVD notrunc;
  
  ProcessTensorRepeat PTR(read_PT);
  PTR.print_info();

  if(PTR.get_n_tot()<1){
    std::cerr<<"PTR.get_n_tot()="<<PTR.get_n_tot()<<"<1!"<<std::endl;
    exit(1);
  }

  TruncatedSVD trunc_tmp=trunc.get_base();
  trunc_tmp.print_info();std::cout<<std::endl;
  if(start_with_forward){
    PTR.sweep_forward(trunc_tmp, false, sweep_initial, verbosity);
  }
  for(int n=0; n<n_sweep; n++){
    std::cout<<"Sweep backward: "<<n<<"/"<<n_sweep<<std::endl;
    PTR.sweep_backward(trunc_tmp, true, true, verbosity);
    for(int m=0; m<n_forward_sweep; m++){
      std::cout<<"Sweep forward: "<<n<<"/"<<n_sweep<<" "<<m<<"/"<<n_forward_sweep+1<<std::endl;
      PTR.sweep_forward(trunc_tmp, true, (m==0), verbosity);
    }
    std::cout<<"Sweep forward: "<<n<<"/"<<n_sweep<<" "<<n_forward_sweep<<"/"<<n_forward_sweep+1<<std::endl;
    PTR.sweep_forward(trunc_tmp, false, (n_forward_sweep<1), verbosity);
  }
  std::cout<<"Final sweep backward"<<std::endl;
  PTR.sweep_backward(trunc_tmp, flip_last, sweep_initial, verbosity);
 
 }catch (DummyException &e){
  return 1;
 }
#ifdef EIGEN_USE_MKL_ALL
  mkl_free_buffers();
#endif

  return 0;
}
