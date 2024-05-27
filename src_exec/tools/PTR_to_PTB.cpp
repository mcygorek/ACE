#include "ProcessTensorBuffer.hpp"
#include "ProcessTensorRepeat.hpp"
#include "TempFileName.hpp"
#include "Parameters.hpp"
#include "GenericSimulation.hpp"
#include "ProcessTensorForwardList.hpp"
#include "DummyException.hpp"

using namespace ACE;
int main(int args, char** argv){

  Parameters param(args, argv, true);

  double dict_zero = param.get_as_double("dict_zero",-1.);
  double dt = param.get_as_double_check("dt");
  double te = param.get_as_double_check("te");
  TimeGrid tgrid(param);
  int n_tot = tgrid.n_tot;
  if(n_tot<1){
    std::cerr<<"n_tot<1!"<<std::endl;
    exit(1);
  }
  int buffer_blocksize = param.get_as_int("buffer_blocksize",-1);
  
  double forward_threshold = param.get_as_double("forward_threshold", 0);
  double backward_threshold = param.get_as_double("backward_threshold", 0);
  
  TruncatedSVD trunc_fw, trunc_bw;
  { Parameters param_fw; 
    param_fw.add_from_prefix("forward",param);
    trunc_fw.setup(param_fw);}
  { Parameters param_bw; 
    param_bw.add_from_prefix("backward",param);
    trunc_bw.setup(param_bw);}


  std::string read_PT = param.get_as_string_check("read_PT");
  std::string write_PT = param.get_as_string_check("write_PT");

  ProcessTensorRepeat PTR(read_PT);
  PTR.set_read_only(true);
  if(PTR.initial.size()<1){
    std::cerr<<"PTR.initial.size()<1!"<<std::endl;
    exit(1);
  }

  ProcessTensorBuffer PTB; 
  PTB.set_new_file(write_PT, buffer_blocksize);
  PTB.create(n_tot);

  PassOn pass_on(PTR.get(0).M.dim_d1);

  for(int n=0; n<n_tot; n++){
    PTB.get(n)=PTR.get(n);
    if(trunc_fw.do_compress()){
      PTB.get(n).sweep_forward(trunc_fw, pass_on, (n==n_tot-1));
    }
  }
  PTB.get(n_tot-1).close_off();

  if(trunc_bw.do_compress()){
    PTB.sweep_backward(trunc_bw,1);
  }

  return 0;
}
