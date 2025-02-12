#ifndef ACE_TRAFOTREE_DEFINED_H
#define ACE_TRAFOTREE_DEFINED_H

#include <Eigen/Dense>
#include "SelectIndices.hpp"
#include "DummyException.hpp"
#include <iostream>
#include <memory>

namespace ACE{
/**
This structure is meant to store local "lossy compression matrices" T and 
their pseudo-inverses T^{-1}. These formally map the full environment 
Liouville space to its most relevant degrees of freedom, i.e, 
Q = T e^{L_E dt} T^{-1}.
In the paper on the inner bonds of PT-MPOs, we showed
that sequential ACE provides these in the form of an half-open MPO. More 
generally (e.g., for the tree-like contraction scheme for ACE), the overall
transformation matrix, which reduces from the full environment Liouville space
to the final inner bonds, has the structure of a binary tree, where the leaves
are local compression matrices.

Thus, we have to implement a binary tree as well as functions to
(a) expand the tree whenever we combine modes (and initialize with Id)
(b) modify the local transformation matrices when sweeping 
(c) storing/reading/writing 
(d) in other classes: do something useful with them (e.g., projecting e^{L_E})  

Preselection might be an issue, as we want to keep the sparsity. In fact, we 
may use the SelectIndices structure to represent the formation of double 
indices. 
This entails some overhead, but we expect this part not to be time critical.
**/

struct TrafoTree{
  std::shared_ptr<TrafoTree> first;
  std::shared_ptr<TrafoTree> second;
  Eigen::MatrixXcd T;
  Eigen::MatrixXcd Tinv;
  SelectIndices S;

  inline bool has_children()const{
    return ((bool)first && (bool)second);
  }
  inline bool is_consistent_T()const{
    return (T.rows()==Tinv.cols() && T.cols()==Tinv.rows());
  }
  inline bool is_consistent_TS()const{
    return (T.cols()==S.size())&&(Tinv.rows()==S.size())&&(T.rows()==Tinv.cols());
  }
  inline int get_dim()const{
    return T.rows();
  }
  inline bool is_consistent_with_children()const{
    if((!(bool)first) && (!(bool)second))return true;
    if((!(bool)first) || (!(bool)second))return false;
    if(!is_consistent_TS())return false;
    return (S.min_dim_first()<=first->get_dim() && 
            S.min_dim_second()<=second->get_dim()); 
  }

  //set T,Tinv to identity matrices with dimension dim
  inline void initialize_T(int dim){
    T=Eigen::MatrixXcd::Identity(dim,dim);
    Tinv=Eigen::MatrixXcd::Identity(dim,dim);
  }
  
  //move contents of this node (including links to first and second) to 
  //the new first child and create second child initialized with identity.
  //void expand(int other_dim);
  
  //move contents of this node (including links to children) to the new first
  //child. Link to "other" as second child.
  //WARNING: "other" will be released (moved)
  void combine_select(std::unique_ptr<TrafoTree> other, 
                      const SelectIndices &select);
  void combine(std::unique_ptr<TrafoTree> other);

  void read(const std::string &filename);
  void write(const std::string &filename)const;

  TrafoTree(int dim=0){
    initialize_T(dim);
  } 
  ~TrafoTree(){
  }
};



}//namespace
#endif
