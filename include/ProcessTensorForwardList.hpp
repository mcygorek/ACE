#ifndef ACE_PROCESS_TENSOR_FORWARD_LIST_DEFINED_H
#define ACE_PROCESS_TENSOR_FORWARD_LIST_DEFINED_H

#include "ProcessTensorForward.hpp"
#include "Eigen_fwd.hpp"
#include "TruncationLayout.hpp"
#include <memory>
#include "Parameters.hpp"
#include "Which_Env_Ops.hpp"
#include "ReadPT_struct.hpp"
#include "ModePropagatorGenerator.hpp"
#include <cstdio>

namespace ACE{

class ProcessTensorForwardList{
public:

  std::vector<std::shared_ptr<ProcessTensorForward> > list;
  std::vector<std::string>   temp_file_list;
  std::vector<ReadPT_struct> temp_expand;
//  std::vector<OuterIndexExtender> extender;

  void complain_if_null()const;
  inline size_t size()const{ return list.size(); }  
  inline void clear(){ list.clear(); temp_expand.clear();}
  void print_info()const;

  void reset();
  bool done()const;
  void load_next();
  std::vector<const ProcessTensorElement *> current_list();

  Eigen::VectorXcd get_rho_reduced(const Eigen::MatrixXcd & state);

  std::vector<std::complex<double> > get_env_reduced(
    const Eigen::MatrixXcd & state, const Which_Env_Ops_List & which_env_ops);


  static std::shared_ptr<ProcessTensorForward> PTptr_from_file(const std::string &str, bool read_only);

  void add_PT(const ReadPT_struct &expand);
  inline void add_PT(const std::string &fname, int front, int back){
    add_PT(ReadPT_struct(fname, front, back));
  }
  void add_PT(Parameters &param);

  void setup2(Parameters &param, 
      std::vector<std::shared_ptr<ModePropagatorGenerator> > & initial_mpgs,
      int setdim=-1, bool print_timings=true);

  inline void setup(Parameters &param, int setdim=-1, bool print_timings=true){
    std::vector<std::shared_ptr<ModePropagatorGenerator> > initial_mpgs;
    setup2(param, initial_mpgs, setdim, print_timings);
  }


  void read(const std::string &fname);
  
  void propagate(Eigen::MatrixXcd & state, bool reverse_order); 
//  void propagate_select(Eigen::MatrixXcd & state, const TruncatedSVD &trunc);

  ProcessTensorForwardList(){}
  ProcessTensorForwardList(const std::string &fname){
    read(fname);
  }
  ProcessTensorForwardList(Parameters &param, int setdim=-1){
    setup(param, setdim);
  }
  ProcessTensorForwardList(Parameters &param, 
      std::vector<std::shared_ptr<ModePropagatorGenerator> > & initial_mpgs){
    setup2(param, initial_mpgs, -1, true);
  }
  ~ProcessTensorForwardList(){
//NOTE: Order of desctruction is crucial!
//    std::cout<<"ProcessTensorForwardList: Destructor called."<<std::endl;
    std::vector<std::shared_ptr<ProcessTensorForward> >().swap(list);
    for(const std::string &str : temp_file_list)std::remove(str.c_str());
  }
};
}//namespace
#endif
