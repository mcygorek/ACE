#ifndef ACE_PROCESS_TENSOR_REPEAT_DEFINED
#define ACE_PROCESS_TENSOR_REPEAT_DEFINED

#include "ProcessTensorForward.hpp"
#include "ProcessTensorElement.hpp"
#include "ProcessTensorBuffer.hpp"
#include "BufferedContainer.hpp"
#include "DiagBB.hpp"
#include "TimeGrid.hpp"
#include "TruncationLayout.hpp"
#include "PreloadHint.hpp"
#include <vector>
#include <climits>


namespace ACE{
/* A process tensor with repeated blocks */

class ProcessTensorRepeat: public ProcessTensorForward{
public: 
  //inherited: int n, n_tot;
  BufferedContainer<ProcessTensorElement> initial;
  BufferedContainer<ProcessTensorElement> repeated;

  virtual int get_n_tot()const{ return INT_MAX; }
//  virtual bool done()const{ return n >= (int)(initial.size()+repeated.size()); }
  virtual bool done()const{ return false; }
  ProcessTensorElement & get(int n_, PreloadHint hint=NoPreload);
  const ProcessTensorElement & get_ro(int n_, PreloadHint hint=NoPreload);

  virtual int get_N_system();
  virtual const ProcessTensorElement * current(){return &get_ro(n, ForwardPreload);}

  static int get_n_mem(const TimeGrid &tgrid);
  static bool fall_back_to_DiagBB_log(const TimeGrid &tgrid, int verbosity=0);

  bool is_read_only(bool val)const;
  void set_read_only(bool val);

  void set_specs(const std::string &write_PT, int blocksize);

  void calculate(DiagBB &diagBB, const TimeGrid &tgrid,
       TruncationLayout trunc, double dict_zero, int verbosity, bool use_log);
 
  void sweep_backward(TruncatedSVD &trunc, bool flip_last, bool sweep_initial, int verbosity=0);
  void sweep_forward(TruncatedSVD &trunc, bool flip_last, bool sweep_initial, int verbosity=0);
  
//  void sweep_final_start_forward(const TruncationLayout &trunc, int verbosity);
  void sweep_final_start_backward(const TruncationLayout &trunc, int verbosity);


  static bool can_read(const std::string &fname_h);
  void read(const std::string &fname_h, bool ro=false);

  virtual void dict_expand(const ReadPT_struct &readPT);

  void set_trivial(int Nsys);

  void print_info(std::ostream &os=std::cout)const;

  ProcessTensorRepeat(int Nsys){
    set_trivial(Nsys);
  }
  ProcessTensorRepeat(const std::string &fname_h="", bool ro=false){
    initial.magicString="PTRI";
    repeated.magicString="PTRR";
    if(fname_h!=""){
      read(fname_h, ro);
    }
  }
  ~ProcessTensorRepeat(){}
};
}
#endif
