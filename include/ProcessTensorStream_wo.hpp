#ifndef ACE_PROCESS_TENSOR_STREAM_WO_DEFINED_H
#define ACE_PROCESS_TENSOR_STREAM_WO_DEFINED_H

#include <fstream>
#include <vector>
#include "ProcessTensorStream_ro.hpp"
#include "ProcessTensorElement.hpp"

namespace ACE{

class ProcessTensorStream_wo{
private: 
  std::ofstream ofs;
  int sz;
  int counter;

public:
  virtual int get_length()const{return sz;}
  virtual int get_counter()const{return counter;}

  void initialize(const std::string &fname, int sz, bool reverse);
  void put(const ProcessTensorElement &e);

  ProcessTensorStream_wo(){}
  ProcessTensorStream_wo(const std::string &fname, int sz_, bool reverse){
    initialize(fname, sz_, reverse);
  }

  void complain_destruct_incomplete()const;
  ~ProcessTensorStream_wo();
};
}//namespace
#endif
