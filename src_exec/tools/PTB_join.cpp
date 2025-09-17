#include "ProcessTensorBuffer.hpp"
#include "TempFileName.hpp"
#include "Parameters.hpp"
#include "GenericSimulation.hpp"
#include "ProcessTensorForwardList.hpp"
#include "DummyException.hpp"
#include "otimes.hpp"

using namespace ACE;
int main(int args, char** argv){

  Parameters param(args, argv, true);

  std::string read_PT = param.get_as_string_check("initial_PT");
//  param.complain_if_not_specified("multi_PT");
  std::vector<std::string> multi_PT = param.get_all_strings("add_PT");
  if(multi_PT.size()<1){
    std::cerr<<"Set at least one 'add_PT'!"<<std::endl;
    exit(1);
  }
 
  double dict_zero = param.get_as_double("dict_zero",-1.);
  int buffer_blocksize = param.get_as_int("buffer_blocksize",-1);
  
  int sweep_n = param.get_as_int("sweep_n",param.get_as_int("intermediate_sweep_n",0));
  int final_sweep_n = param.get_as_int("final_sweep_n",0);
  int verbosity = param.get_as_int("verbosity", 2); 
 
  bool only_join = param.get_as_bool("only_join", false);
  bool use_select = param.get_as_bool("use_select", true);
  bool force_NF = param.get_as_bool("force_NF", false);

  TruncationLayout trunc(param);
  TruncatedSVD notrunc;
  
/*
  std::string write_PT = param.get_as_string_check("write_PT");
  ProcessTensorBuffer PTB;
  PTB.set_new_file(write_PT, buffer_blocksize);
  {ReadPT_struct read_struct(param, "read_PT", "read_PT_expand");
  std::cout<<"Reading: '"<<read_struct.fname<<"'"<<std::endl;
  PTB.copy_content(read_struct);
  }
*/
  ProcessTensorBuffer PTB(read_PT);

  if(PTB.get_n_tot()<1){
    std::cerr<<"PTB.get_n_tot()="<<PTB.get_n_tot()<<"<1!"<<std::endl;
    exit(1);
  }
  for(int i=0; i<(int)multi_PT.size(); i++){
    const std::string & mPT = multi_PT[i];
    ProcessTensorBuffer PTB1(mPT,true);
    if(PTB1.get_n_tot()<1){
      std::cerr<<"'"<<mPT<<"': PTB1.get_n_tot()="<<PTB1.get_n_tot()<<"<1!"<<std::endl;
      throw DummyException();
    }
    if(!PTB1.get(0).is_forwardNF()){
      std::cerr<<"'"<<mPT<<"': not in forward normal form. Please run 'PTB_sweep_forward' before running PTB_join!"<<std::endl;
      throw DummyException();
    }
  }


  for(int i=0; i<(int)multi_PT.size(); i++){
    const std::string & mPT = multi_PT[i];
    ProcessTensorBuffer PTB1(mPT,true);

//    TruncatedSVD trunc_first=trunc.get_forward(i,multi_PT.size());
    TruncatedSVD trunc_first;
    if(!PTB.get(0).is_forwardNF() || force_NF){
      std::cout<<"Sweeping forward first PT"<<std::endl; // '"<<write_PT<<"'"<<std::endl;
      trunc_first.print_info(); std::cout<<std::endl;
      PTB.sweep_forward(trunc_first, verbosity);
    }

    std::cout<<"Joining and sweeping forward"<<std::endl;
    TruncatedSVD trunc_select=trunc.get_backward(i,multi_PT.size(),true);
    TruncatedSVD trunc_bw=trunc.get_backward(i,multi_PT.size());
    trunc_select.print_info(); std::cout<<std::endl;

    if(only_join){
      for(int n=PTB.get_n_tot()-1; n>=0; n--){
        std::cout<<"n="<<n<<"/"<<PTB.get_n_tot()<<std::endl;
        if(n>0&&n==param.get_as_int("n_break"))break;

        ProcessTensorElement & e = PTB.get(n, BackwardPreload);
        ProcessTensorElement & e2 = PTB1.get(n, BackwardPreload);
        std::cout<<"e.forwardNF="<<e.forwardNF.transpose()<<std::endl;
        std::cout<<"e2.forwardNF="<<e2.forwardNF.transpose()<<std::endl;

/*        std::cout<<"e.M:"<<std::endl;
        for(int i=0; i<e.M.dim_i; i++){
          std::cout<<"i="<<i<<std::endl;
          std::cout<<e.M.get_Matrix_d1_d2(i)<<std::endl;
        }
        std::cout<<"e2.M:"<<std::endl;
        for(int i=0; i<e2.M.dim_i; i++){
          std::cout<<"i="<<i<<std::endl;
          std::cout<<e2.M.get_Matrix_d1_d2(i)<<std::endl;
        }*/
//        e.join_thisfirst(e2);
        e.accessor.join_thisfirst(e2.accessor, e.M, e2.M);
        e.closure=Vector_otimes(e.closure, e2.closure);
        e.env_ops.join(e2.env_ops);
        e.clearNF();
/*
        std::cout<<"e.M:"<<std::endl;
        for(int i=0; i<e.M.dim_i; i++){
          std::cout<<"i="<<i<<std::endl;
          std::cout<<e.M.get_Matrix_d1_d2(i)<<std::endl;
        }
*/
        std::cout<<std::endl;
      }
    }else if(use_select){
      PTB.join_select_and_sweep_backward(PTB1, trunc_select, trunc_bw, verbosity);
    }else{ 
      PTB.join_and_sweep_backward(PTB1, trunc_bw, verbosity);
    }


    PTB.sweep_intermediate_or_final_start_forward(trunc, i, multi_PT.size(), verbosity);
 }


 

#ifdef EIGEN_USE_MKL_ALL
  mkl_free_buffers();
#endif

  return 0;
}
