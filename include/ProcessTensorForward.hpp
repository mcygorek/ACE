#ifndef ACE_PROCESS_TENSOR_FORWARD_DEFINED_H
#define ACE_PROCESS_TENSOR_FORWARD_DEFINED_H
#include "ProcessTensorElement.hpp"

namespace ACE{

/* Abstraction: sequential access to PT elements (in forward direction)
*/
class ProcessTensorForward{
public:
  int n;
  int n_tot;

  virtual int get_n_tot()const{ return n_tot; }
  virtual void reset(){ n=0; }
  virtual bool done()const{ return n>=get_n_tot(); }
  virtual void load_next(){ n++; }
  virtual void load(int l){ 
    if(l==n){return;}
    if(l>n){ for(int i=n; i<l; i++){load_next();} }
    reset(); for(int i=0; i<l; i++){load_next();}
  }

  virtual void dict_expand(const ReadPT_struct &readPT) = 0;
  virtual int get_N_system() = 0;
  virtual const ProcessTensorElement * current() = 0;
  ProcessTensorForward(){
    n=0;
    n_tot=0; 
  }
};
}//namespace
#endif
