#ifndef ACE_PROCESS_TENSOR_DEFINED_H
#define ACE_PROCESS_TENSOR_DEFINED_H

#include "ProcessTensorElement.hpp"
#include <vector>
#include "TimeGrid.hpp"
#include "ModePropagatorGenerator.hpp"

namespace ACE{

enum COMBINE_MODE{ mode_sequential, mode_select };

/*
class ProcessTensor_abstract{
public:
  virtual const ProcessTensorElement & current()const=0;
  virtual ProcessTensorElement & current()=0;
};

class ProcessTensor_read_forward: virtual public ProcessTensor_abstract{
public:
  virtual bool is_last()const=0;
  virtual void move_forward()=0;
};
class ProcessTensor_read_backward: virtual public ProcessTensor_abstract{
public:
  virtual bool is_first()const=0;
  virtual void move_backward()=0;
};
*/

class ProcessTensor{
public:

  std::vector<ProcessTensorElement> elements;

  inline size_t size()const{
    return elements.size();
  }
  inline void resize(int nmax){
    std::vector<ProcessTensorElement> empty(nmax);
    elements.swap(empty);
  }
  inline const ProcessTensorElement & back()const{return elements.back();}
  inline ProcessTensorElement & back(){return elements.back();}

  virtual const ProcessTensorElement & operator[](size_t i)const;
  virtual ProcessTensorElement & operator[](size_t i);

  int get_max_dim()const;
  void check_compatible(const ProcessTensor &other, int factor=1);

  void sweep_forward(const TruncatedSVD &trunc);
  void sweep_backward(const TruncatedSVD &trunc);
  void join_halfdt(ProcessTensor &other);
  void join_and_sweep_halfdt(ProcessTensor &other, const TruncatedSVD &trunc);
//  void join_and_sweep_select(ProcessTensor &other, const TruncatedSVD &trunc, int verbosity);

  //requires both PTs in forward normal form
  void join_select_and_sweep_backward(ProcessTensor &other, const TruncatedSVD &trunc);

 
  void set_from_ModePropagator(ModePropagator &mprop, const TimeGrid &tgrid, double dict_zero=0);
  void set_trivial(int Nsys, const TimeGrid &tgrid);
  void calculate_closures();
  
  void print_dims(std::ostream &ofs=std::cout)const;

  void add_modes(COMBINE_MODE combine_mode, 
          ModePropagatorGenerator &mpg, const TimeGrid &tgrid, 
          const TruncatedSVD & trunc, int intermediate_sweep_n,
          double dict_zero, int verbosity);


  static std::pair<int, bool> read_header(std::ifstream &is, const std::string & context="");
  static void write_header(std::ofstream &os, int sz, bool reverse);
  static int read_size(const std::string &fname);

  void read_binary(std::ifstream &is, const std::string &context="");
  void read_binary(const std::string &fname);
  //reverse: If set, the elements are stored in reverse order (used for line sweeps with storage on hard disk):
  void write_binary(std::ofstream &os, bool reverse=false)const;
  void write_binary(const std::string &fname, bool reverse=false)const;


  inline ProcessTensor(){}
  inline ProcessTensor(int Nsys, const TimeGrid &tgrid){
    set_trivial(Nsys, tgrid);
  }
  inline ProcessTensor(const std::string& fname){
    read_binary(fname);
  }
  inline ProcessTensor(ModePropagator &mprop, const TimeGrid &tgrid, double dict_zero=0){ 
    set_from_ModePropagator(mprop, tgrid, dict_zero);
  }
};


}//namespace
#endif
