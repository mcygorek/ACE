#include "ProcessTensorBuffer.hpp"
#include "Parameters.hpp"
#include "TruncatedSVD.hpp"
#include "AddPT.hpp"
#include <iostream>

using namespace ACE;

int main(int args, char **argv){

 try{

  Parameters param(args, argv, true);
  TruncatedSVD trunc(param);

  std::string read_PT=param.get_as_string_check("read_PT");
//  std::string write_PT=param.get_as_string_check("write_PT");
//  int buffer_blocksize=param.get_as_int("buffer_blocksize",-1);

  ProcessTensorBuffer PTB(read_PT);
  ProcessTensorBuffer PTB2;
//  PTB2.set_new_file(write_PT, buffer_blocksize);
//  PTB2.set_new_temporary(PTB);

  bool ret=PTB.split_inner(PTB2, PTB.get_n_tot()/2);
//  bool ret=PTB.split_inner_and_sweep_fbf(PTB2, trunc, PTB.get_n_tot()/2, 1);
  if(ret){
    std::cout<<"Split successful."<<std::endl;
  }else{
    std::cout<<"Split unsuccessful."<<std::endl;
  }

  PTB.additive_join(PTB2);
//  PTB.calculate_closures();

  std::cout<<"Done."<<std::endl;

 }catch (DummyException &e){
 }

  return 0;
}
