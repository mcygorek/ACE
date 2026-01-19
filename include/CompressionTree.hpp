#ifndef ACE_COMPRESSIONTREE_DEFINED_H
#define ACE_COMPRESSIONTREE_DEFINED_H

#include <Eigen/Dense>
#include "SelectIndices.hpp"
#include "DummyException.hpp"
#include "BinaryReader.hpp"
#include "PassOn.hpp"
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

This only works for the method ACE with or without use_combine_tree!

There should be one global share_ptr of the full tree, and the 
ProcessTensorElementAccessor should contain a share_ptr to it if the cut 
is to the left or to the right of the corresponding PT element.
If the PT element finds itself to be combined or compressed, it should
call back to functions in the tree structure. 

This entails a problem when copying PT elements: Either one has to deep 
copy the full tree or one has to move to tree to only one of the elements.
Here, we choose not to copy the tree automatically but leave this to
the join/sweep functions in ProcessTensorElementAccessor. 

Also, the pointer has to be shared_ptr, because elements left and right of
the cut have to be able to point to the same CompressionTree. This means, one should
ABSOLUTELY NOT copy ProcessTensorElements when also using CompressionTrees, because
two elements point to the same tree and most likely cause inconsistencies.

Concretely, we have to implement here:
(a) functions and elements relating to the tree structure 
(a) expand the tree whenever we combine modes (and initialize with Id)
(b) modify the local transformation matrices when sweeping 
(c) storing/reading/writing

Doing something with it, such as projecting e^{L_E}, is left for other classes
that take the complete tree as input.

Preselection might be an issue, as we want to keep the sparsity. In fact, we 
may use the SelectIndices structure to represent the formation of double 
indices. 
This entails some overhead, but we expect this part not to be time critical.
**/

struct CompressionTree{
  std::shared_ptr<CompressionTree> first;
  std::shared_ptr<CompressionTree> second;
  Eigen::MatrixXcd T;
  SelectIndices S;

  inline bool has_children()const{
    return ((bool)first && (bool)second);
  }
  inline bool is_consistent_TS()const{
    return (T.cols()==S.size());
  }
  inline int get_dim()const{
    return T.rows();
  }
  inline bool is_consistent_with_children()const{
    if(!has_children())return true;
    if(!is_consistent_TS())return false;
    return (S.min_dim_first()<=first->get_dim() && 
            S.min_dim_second()<=second->get_dim()); 
  }

  //set T to identity matrices with dimension dim
  inline void initialize_T(int dim){
    T=Eigen::MatrixXcd::Identity(dim,dim);
  }

  inline void print_info(std::ostream &out=std::cout)const{
    out<<"TTree: T.rows()="<<T.rows()<<" T.cols()="<<T.cols()<<" S.size()="<<S.size(); 
    if(first)out<<" first=true";else out<<" first=false";
    if(second)out<<" second=true";else out<<" second=false";
    out<<std::endl;
  }    
  inline void print_info_traverse(const std::string &pad="", std::ostream &out=std::cout)const{
    out<<pad<<"T: "<<T.rows()<<","<<T.cols()<<" S: "<<S.size()<<","<<S.min_dim_first()<<","<<S.min_dim_second()<<std::endl;
    if(first){out<<"first:"<<std::endl;
              first->print_info_traverse(pad+".", out);
    }else{out<<"first: none."<<std::endl;}
    if(second){out<<"second:"<<std::endl;
               second->print_info_traverse(pad+".", out);
    }else{out<<"second: none."<<std::endl;}
  }
  inline void print_info_file(const std::string &fname)const{
    std::ofstream ofs(fname);
    print_info_traverse("", ofs);
  }


  //move contents of this node (including links to first and second) to 
  //the new first child and create second child initialized with identity.
  //void expand(int other_dim);
  
  //move contents of this node (including links to children) to the new first
  //child. Link to "other" as second child.
  //WARNING: "other" will be released (moved)
  void combine_select(std::shared_ptr<CompressionTree> & other, 
                      const SelectIndices &select);
  void combine(std::shared_ptr<CompressionTree> & other);


  void read(std::istream &ifs, const std::string &context="");
  void write(std::ostream &ofs)const;

  inline void read(const std::string &filename){
    std::unique_ptr<std::ifstream> ifs = open_file_check(filename);
    read(*(ifs.get()),filename);
  }
  inline void write(const std::string &filename)const{
    std::ofstream ofs(filename);
    write(ofs);
  }

  CompressionTree(int dim=0){
    initialize_T(dim);
  } 
  CompressionTree(std::istream &ifs, const std::string &context=""){
    read(ifs, context);
  }
  CompressionTree(const std::string &filename){
    read(filename);
  }
  ~CompressionTree(){
  }
};


void sweep_forward_pre(int n, int TTree_at, \
                       std::shared_ptr<CompressionTree> & TTree, \
                       std::shared_ptr<CompressionTree> & TTree_inv, \
                       const PassOn &pass_on);

void sweep_forward_post(int n, int TTree_at, \
                       std::shared_ptr<CompressionTree> & TTree, \
                       std::shared_ptr<CompressionTree> & TTree_inv, \
                       const PassOn &pass_on);

void sweep_backward_pre(int n, int TTree_at, \
                       std::shared_ptr<CompressionTree> & TTree, \
                       std::shared_ptr<CompressionTree> & TTree_inv, \
                       const PassOn &pass_on);

void sweep_backward_post(int n, int TTree_at, \
                       std::shared_ptr<CompressionTree> & TTree, \
                       std::shared_ptr<CompressionTree> & TTree_inv, \
                       const PassOn &pass_on);
}//namespace
#endif
