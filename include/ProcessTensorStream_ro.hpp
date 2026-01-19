#ifndef ACE_PROCESS_TENSOR_STREAM_RO_DEFINED_H
#define ACE_PROCESS_TENSOR_STREAM_RO_DEFINED_H

#include <fstream>
#include <vector>
#include "ProcessTensorElement.hpp"
#include "ProcessTensorForward.hpp"
#include <memory>

namespace ACE{

class ProcessTensorStream_ro: public ProcessTensorForward{
private:
  std::ifstream ifs;
  std::vector<std::streampos> s_pos;

  std::shared_ptr<ProcessTensorElement> fwd_buffer;
public:
  inline size_t size()const{return s_pos.size();}
  void complain_if_empty()const;

  ProcessTensorElement get(int n);

  inline ProcessTensorElement get_first(){
    complain_if_empty();
    return get(0);
  }
  inline ProcessTensorElement get_last(){
    complain_if_empty();
    return get((int)size()-1);
  }
  inline ProcessTensorElement get_as_last(int n){
    complain_if_empty();
    ProcessTensorElement e=get(n);
    e.close_off();
    return e;
  }

  virtual int get_n_tot()const;
  virtual void dict_expand(const ReadPT_struct &readPT);
  virtual int get_N_system();
  virtual const ProcessTensorElement * current();
  
  void process_file(const std::string &fname, bool complain=false);

  ProcessTensorStream_ro(){}
  ProcessTensorStream_ro(const std::string &fname, bool complain=false){
    process_file(fname, complain);
  }
};

}//namespace
#endif
