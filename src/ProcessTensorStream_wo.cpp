#include "ProcessTensorStream_wo.hpp"
#include "ProcessTensor.hpp"
#include <iostream>
#include "DummyException.hpp"


namespace ACE{

void ProcessTensorStream_wo::initialize(const std::string &fname, int sz_, bool reverse){
  sz=sz_;
  counter=0;
  if(ofs.is_open())ofs.close();
  ofs.clear();
  ofs.open(fname.c_str());
  ProcessTensor::write_header(ofs, sz, reverse);
}  
void ProcessTensorStream_wo::put(const ProcessTensorElement &e){
  if(sz==0 || !ofs.is_open() || !ofs.good()){
    std::cerr<<"ProcessTensorStream_wo::put: sz==0 || !ofs.is_open() || !ofs.good()!"<<std::endl;
    throw DummyException();
  }
  if(counter>=sz){
    std::cerr<<"ProcessTensorStream_wo::put: counter>=sz ("<<counter<<" vs. "<<sz<<")!"<<std::endl;
    throw DummyException();
  }
  e.write_binary(ofs);
  counter++;
}

void ProcessTensorStream_wo::complain_destruct_incomplete()const{
  if(ofs.is_open()){
    if(counter!=sz){
      std::cerr<<"Destruction of incomplete ProcessTensorStream_wo with counter!=sz-1 ("<<counter<<"/"<<sz<<")!"<<std::endl;
      throw DummyException();
    }
  }
}
ProcessTensorStream_wo::~ProcessTensorStream_wo(){
  complain_destruct_incomplete();
}




}//namespace
 
