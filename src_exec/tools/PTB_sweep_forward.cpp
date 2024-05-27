#include "ProcessTensorBuffer.hpp"
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
  
  int sweep_n = param.get_as_int("sweep_n",param.get_as_int("intermediate_sweep_n",0));
  int verbosity = param.get_as_int("verbosity", 2); 

  TruncationLayout trunc(param);
  trunc.print_info(); std::cout<<std::endl;
  TruncatedSVD notrunc;
  
  ProcessTensorBuffer PTB(read_PT);

  if(PTB.get_n_tot()<1){
    std::cerr<<"PTB.get_n_tot()="<<PTB.get_n_tot()<<"<1!"<<std::endl;
    exit(1);
  }

  for(int loop=0; loop<sweep_n; loop++){
    TruncatedSVD trunc_tmp=trunc.get_base();
    std::cout<<"loop="<<loop<<"/"<<sweep_n<<std::endl;
    std::cout<<"Sweeping forward"<<std::endl;
    trunc_tmp.print_info();std::cout<<std::endl;
    PTB.sweep_forward(trunc_tmp, verbosity);

    trunc_tmp=trunc.get_base();
    std::cout<<"Sweeping backward"<<std::endl;
    trunc_tmp.print_info();std::cout<<std::endl;
    PTB.sweep_backward(trunc_tmp, verbosity);
  }

  TruncatedSVD trunc_tmp=trunc.get_base();
  std::cout<<"Sweeping forward"<<std::endl;
  trunc_tmp.print_info();std::cout<<std::endl;
  PTB.sweep_forward(trunc_tmp, verbosity);
 
 }catch (DummyException &e){
  return 1;
 }
#ifdef EIGEN_USE_MKL_ALL
  mkl_free_buffers();
#endif

  return 0;
}
