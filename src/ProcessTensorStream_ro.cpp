#include "ProcessTensorStream_ro.hpp"
#include "ProcessTensor.hpp"
#include "DummyException.hpp"

namespace ACE{

void ProcessTensorStream_ro::complain_if_empty()const{
  if(size()<1){
    std::cerr<<"ProcessTensorStream_ro: PT is empty!"<<std::endl;
    throw DummyException();
  }
}

ProcessTensorElement ProcessTensorStream_ro::get(int n){
  if(n<0||n>=(int)s_pos.size()){
    std::cerr<<"ProcessTensorStream_ro::get: Out of bounds: "<<n<<"/"<<s_pos.size()<<std::endl;
    throw DummyException();
  }
  if(!ifs.seekg(s_pos[n]).good()){
    std::cerr<<"ProcessTensorStream_ro::get("<<n<<"): seekg failed!"<<std::endl;
    throw DummyException();
  }
  ProcessTensorElement e;
  e.read_binary(ifs);
  if(!ifs.seekg(s_pos[n]).good()){
    std::cerr<<"ProcessTensorStream_ro::get("<<n<<"): read_binary failed!"<<std::endl;
    throw DummyException();
  }
  return e;
}

int ProcessTensorStream_ro::get_n_tot()const{
  return s_pos.size();
}

int ProcessTensorStream_ro::get_N_system(){
  int this_n=ProcessTensorForward::n; 
  if(this_n>=n_tot){ this_n=0; }
  if(n_tot<1){
    std::cerr<<"ProcessTensorStream_ro::get_N_system(): n_tot<1!"<<std::endl;
    throw DummyException();
  }
  const ProcessTensorElement &e=get(this_n);
  return e.get_N();
}
const ProcessTensorElement * ProcessTensorStream_ro::current(){
  fwd_buffer = std::shared_ptr<ProcessTensorElement>(new ProcessTensorElement(get(ProcessTensorForward::n)));
  return fwd_buffer.get();
}

void ProcessTensorStream_ro::dict_expand(const ReadPT_struct &readPT){
  std::cerr<<"Cannot expand system dimension of ProcessTensorStream_ro!"<<std::endl;
  readPT.print_info(std::cerr);
  throw DummyException();
}

void ProcessTensorStream_ro::process_file(const std::string &fname, bool complain){
  if(ifs.is_open())ifs.close();
  ifs.clear();
  ifs.open(fname.c_str());

  std::pair<int, bool> header = ProcessTensor::read_header(ifs, fname);
  int sz=header.first;
  bool reverse=header.second;

  s_pos.resize(sz);
  ProcessTensorElement dummy;
  if(reverse){
    for(int n=(int)s_pos.size()-1; n>=0; n--){
      s_pos[n]=ifs.tellg();
      dummy.read_binary(ifs);
    }
  }else{
    for(int n=0; n<(int)s_pos.size(); n++){
      s_pos[n]=ifs.tellg();
      dummy.read_binary(ifs);
    }
  }

  n=0; n_tot=sz; 
  if(complain)complain_if_empty();
}

}
