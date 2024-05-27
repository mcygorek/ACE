#ifndef ACE_PROCESS_TENSOR_ELEMENT_ACCESSOR_DEFINED_H
#define ACE_PROCESS_TENSOR_ELEMENT_ACCESSOR_DEFINED_H

#include "MPS_Matrix.hpp"
#include "Eigen_fwd.hpp"
#include <vector>
#include "IF_OD_Dictionary.hpp"
#include "SelectIndices.hpp"
#include "HilbertSpaceRotation.hpp"

namespace ACE{

class ProcessTensorElementAccessor{
public:
  typedef std::vector<std::vector<std::pair<int, int> > > VVPI;

  //Note: this may update "dict" and leave "M" inconsistent
  VVPI join_thisfirst_indices(const ProcessTensorElementAccessor &other);  
  VVPI join_thissecond_indices(const ProcessTensorElementAccessor &other);  
  void join(const VVPI &i_list, MPS_Matrix & M, const MPS_Matrix & M2);
  void join_select_indices(const VVPI &i_list,
               MPS_Matrix &M1, const MPS_Matrix &M2,
               const ProcessTensorElementAccessor & acc_other,
               const SelectIndices & k_left, const SelectIndices & k_right);



  IF_OD_Dictionary dict;

  inline void swap(ProcessTensorElementAccessor &other){
    dict.swap(other.dict);
  }
  //returns system dimension
  inline int get_N()const{
    return dict.get_N();
  } 
  //quit if mismatch in system dimension
  void check_N(int dim)const;

  inline void dict_expand(const ReadPT_struct &readPT){
    dict.expand(readPT);
  }

  void set_trivial(int N, MPS_Matrix &M);

  inline void set_from_dict(const IF_OD_Dictionary & dct){
    dict=dct;
  }

  inline void expand_DiagBB(const DiagBB &diagBB){
    dict.expand_DiagBB(diagBB);
  }
  void apply_HilbertSpaceRotation(const HilbertSpaceRotation &hs_rot, MPS_Matrix &M, double dict_zero=0.);

  //combination of system indices describing trace preservation for the process of calculating closures.
  std::pair<double, std::vector<int> > closure_indices()const;

  //join two elements and replace this one. 
  //The outer indices are multiplied as indicated by the accessors where 
  //this element is applied first (the other one second):
  virtual void join_thisfirst(const ProcessTensorElementAccessor &other,
                      MPS_Matrix & M_this, const MPS_Matrix & M_other);

  //Same but this element is applied second (the other one first):
  virtual void join_thissecond(const ProcessTensorElementAccessor &other,
                     MPS_Matrix & M_this, const MPS_Matrix & M_other);

  //Join as in symmetric Trotter combination: d2(left) connected to d1(right)
  void join_symmetric(
         const ProcessTensorElementAccessor & acc_L,
         const ProcessTensorElementAccessor & acc_R,
         MPS_Matrix & M, const MPS_Matrix & M_L, const MPS_Matrix & M_R);


  //as above, however, restrict double indices to certain combinations
  void join_select_indices_alternate(
               int n, MPS_Matrix &M1, const MPS_Matrix &M2,
               const ProcessTensorElementAccessor & acc_other,
               const SelectIndices & k_left, const SelectIndices & k_right);

  //void join_meanfield(IF_OD_Dictionary other_dict, MPS_Matrix & M, MPS_Matrix M2);

  //used for "coarse graining": join along both, outer and inner indices:
//  void join_thisfirst_sameinner(ProcessTensorElement &e1, const ProcessTensorElementAccessor &e2);  //-> moved to ProcessTensorElement


  void propagate(Eigen::MatrixXcd &state, const MPS_Matrix &M, int dim1_front,
                 const ReadPT_struct &expand=ReadPT_struct())const;
  
  void read_binary(std::istream &is);
  void write_binary(std::ostream &os)const;
};

}//namespace


#endif
