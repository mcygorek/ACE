#include "ProcessTensorBufferSpec.hpp"
#include "DummyException.hpp"
#include "TempFileName.hpp"
#include <string>

namespace ACE{

std::string ProcessTensorBufferSpec::get_fname(int n){
  if(use_single_file){
    if(n!=0){
      std::cerr<<"ProcessTensorBufferSpec: get_fname(n) called with use_single_file=true but n="<<n<<"!=0!"<<std::endl;
      throw DummyException(); 
    }
    return fname_header;
  }else{
    return fname_header+"_"+std::to_string(n);
  }
}

}//namespace
