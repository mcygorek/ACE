#ifndef ACE_PROCESS_TENSOR_ELEMENT_DEFINED_H
#define ACE_PROCESS_TENSOR_ELEMENT_DEFINED_H

#include <memory>
//#include "InfluenceFunctional_OD.hpp"
#include "ProcessTensorElementAccessor.hpp"
#include "PassOn.hpp"
#include "EnvironmentOperators.hpp"
#include "MPS_Matrix.hpp"
#include "ModePropagator.hpp"
#include "TruncatedSVD.hpp"
#include "HilbertSpaceRotation.hpp"
#include "DiagBB.hpp"

/*
  high-level object describing an individual process tensor element
*/

namespace ACE{
class InfluenceFunctional_OD;

struct ProcessTensorElement{

  //Metadata:
 // double dt; //timestep

  //handles how to access and combine outer indices, e.g. dict,...
  ProcessTensorElementAccessor accessor;
 
  //closure+env_ops:
  Eigen::VectorXcd closure;
  EnvironmentOperators env_ops;

  //optional data: 
  Eigen::VectorXd forwardNF, backwardNF;

  //Actual MPS matrix storage
  MPS_Matrix M;

  //only relevant during creation: 
  //std::shared_ptr<Trafo_Chain> link_left_TC;
  //std::shared_ptr<Trafo_Chain> link_right_TC;
  //

  inline void swap(ProcessTensorElement &e){
    accessor.swap(e.accessor);
    closure.swap(e.closure);
    env_ops.swap(e.env_ops);
    forwardNF.swap(e.forwardNF); 
    backwardNF.swap(e.backwardNF);
    M.swap(e.M);
  }

  //returns system dimension
  inline int get_N()const{  return accessor.get_N();  }
  //quit if mismatch in system dimension
  inline void check_N(int dim)const{ accessor.check_N(dim); }
  //
  void check_consistency(const std::string & context="")const;

 
  //check if normal form
  bool is_forwardNF()const;
  bool is_backwardNF()const;
  void printNF(std::ostream &ofs=std::cout)const;
  //remove normal form status
  void clearNF();
  

  //apply to concrete state:
  inline void propagate(Eigen::MatrixXcd &state, int dim1_front,
                        const ReadPT_struct &expand=ReadPT_struct())const{
    accessor.propagate(state, M, dim1_front, expand);
  }

  //join two elements and replace this one. 
  //The outer indices are multiplied as indicated by the accessors where 
  //this element is applied first (the other one second):
  void join_thisfirst(const ProcessTensorElement &other); 

  //Same but this element is applied second (the other one first):
  void join_thissecond(const ProcessTensorElement &other); 

  inline void join(const ProcessTensorElement &other, bool this_second){
    if(this_second){
      join_thissecond(other);
    }else{
      join_thisfirst(other);
    }
  } 

  //Join as in symmetric Trotter combination: d2(left) connected to d1(right)
  void join_symmetric(const ProcessTensorElement &other_left, 
                      const ProcessTensorElement &other_right); 
 
  //join along both, outer and inner indices: used, e.g., for coarse graining
  void join_thisfirst_sameinner(const ProcessTensorElement &e2);

 
  //local compression: returned value is the passed-on matrix
  void sweep_forward(const TruncatedSVD &trunc, PassOn &pass_on, bool is_last);
  void sweep_backward(const TruncatedSVD &trunc, PassOn &pass_on, bool is_last);
  
  //canonicalization using QR; so far no truncation implemented
  void sweep_forward_QR(const TruncatedSVD &trunc, PassOn &pass_on, bool is_last);
  void sweep_backward_QR(const TruncatedSVD &trunc, PassOn &pass_on, bool is_last);

  //pair-wise compression to escape local minima
  void sweep_pair_forward(ProcessTensorElement &e2, const TruncatedSVD &trunc);
  void sweep_pair_backward(ProcessTensorElement &e2, const TruncatedSVD &trunc);


  //Block combination:
  SelectIndices get_forwardNF_selected_indices(
      const ProcessTensorElement &other, const TruncatedSVD &trunc)const;

  SelectIndices get_backwardNF_selected_indices(
      const ProcessTensorElement &other, const TruncatedSVD &trunc)const;


  void join_selected(int n, const ProcessTensorElement &other,
       const SelectIndices & k_list_left, const SelectIndices & k_list_right);

  void join_average_selected(const ProcessTensorElement &other,
       const SelectIndices & k_list_left, const SelectIndices & k_list_right);

  //clears everything
  void clear();
  void set_trivial(int N_sys);

  void set_from_ModePropagator(ModePropagator &mprop, double ta, double dt, double dict_zero=0);

  void set_from_InfluenceFunctional_OD(const InfluenceFunctional_OD &IF, int n);

  IF_OD_Dictionary detect_dict(double dict_zero)const;
  void reduce_to_dict(const IF_OD_Dictionary &dict);
  void expand_from_dict();
  void expand_DiagBB(const DiagBB &diagBB);
  void expand_space_front(int N_front);
  void expand_space_back(int N_back);
  void apply_HilbertSpaceRotation(const HilbertSpaceRotation &hs_rot, double dict_zero=0.);

  void calculate_closure(const ProcessTensorElement *last);

  //reduce right dim to 1 by multiplying with closure
  void close_off();
 
  
  inline void print_dims(std::ostream &ofs=std::cout)const{
    return M.print_dims(ofs);
  }
  void print_debug(std::ostream &ofs=std::cout)const;

  void read_binary(std::istream &is);
  void write_binary(std::ostream &os)const;
  

  ProcessTensorElement(){}
  ProcessTensorElement(int N_sys){
    set_trivial(N_sys);
  }
  ProcessTensorElement(const InfluenceFunctional_OD &IF, int n){
   set_from_InfluenceFunctional_OD(IF, n);
  }
};

}//namespace

#endif
