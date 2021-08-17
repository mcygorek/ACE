#ifndef ACE_SPACE_LAYOUT_DEFINED_H
#define ACE_SPACE_LAYOUT_DEFINED_H

#include "CheckMatrix.h"
#include <vector>

/*  Purpose:
Defines 'physical' layout, i.e., meaning of indices on matrices.

For example, the system-environment Hamiltonian is a matrix on Hilbert space
(H_S otimes H_E). The full propagator maps from the Liouville space
((H_S otimes H_E) otimes (H_S otimes H_E)) = (H_S otimes H_E)^2 onto itself.
If we want to combine system-environment Hamiltonians for two different
environments, we have to expand Hilbert spaces, e.g.,
(H_S otimes H_E2) -> (H_S otimes H_E1 otimes H_E2)  or the respective 
Liouville spaces (H_S otimes H_E2)^2 -> (H_S otimes H_E1 otimes H_E2)^2

All of the above is typically represented in terms of matrices, whose
indices are collective indices. Here, we between Kronecker product indices
(single int) and the indices corresponding to different subspaces (int vector).

NOTE: This implementation favours correctness over speed!

*/

class ProductSpace{
private:
  std::vector<int> dims; //dimensions of subspaces
  int total_size; //keep consistent with every change in dims!
public:
  int get_dim(int l)const{
    check_bounds(l, dims.size(), "ProductSpace::get_dim(int l): l");
    return dims[l];
  }
  int get_total_size()const{
    return total_size;
  }
  size_t get_dims_size()const{
    return dims.size();
  }
  void erase(int l){
    check_bounds(l, dims.size(), "ProductSpace::erase(int l): l");
    total_size/=dims[l];
    dims.erase(dims.begin()+l);
  }
 
  void check_validity(const std::string &name=""){
    if(dims.size()<1){
      std::stringstream ss;
      ss<<"ProductSpace "; if(name!=""){ss<<"'"<<name<<"' ";}; 
      ss<<"has no dimension set. Not initialized?"<<std::endl;
      throw(std::runtime_error(ss.str()));
    }
    for(size_t i=0; i<dims.size(); i++){
      if(dims[i]<1){
        std::stringstream ss;
        ss<<"ProductSpace "; if(name!=""){ss<<"'"<<name<<"' ";}; 
        ss<<"dims["<<i<<"]="<<dims[i]<<". Must be at least 1!"<<std::endl;
        throw(std::runtime_error(ss.str()));
      }
    }
  }
  void recalculate_total_size(const std::string &name=""){
    check_validity(name);
    total_size=1;
    for(size_t l=0; l<dims.size(); ++l){
      total_size*=dims[l];
    }
  }
   
  ProductSpace square(const std::string &name="")const{
    std::vector<int> v(dims); 
    v.insert(v.end(), dims.begin(), dims.end());
    return ProductSpace(v,name);
  }
 
  ProductSpace(const std::vector<int> &dims, const std::string &name="") 
   : dims(dims){
    recalculate_total_size(name);
  }
  ProductSpace(int i0, int i1, int i2, int i3, const std::string &name=""){
    dims.resize(4);
    dims[0]=i0;  
    dims[1]=i1;  
    dims[2]=i2;  
    dims[3]=i3;  
    recalculate_total_size(name);
  }
  ProductSpace(int i0, int i1, int i2, const std::string &name=""){
    dims.resize(3);
    dims[0]=i0;  
    dims[1]=i1;  
    dims[2]=i2;  
    recalculate_total_size(name);
  }
  ProductSpace(int i0, int i1, const std::string &name=""){
    dims.resize(2);
    dims[0]=i0;  
    dims[1]=i1;  
    recalculate_total_size(name);
  }
  ProductSpace(int i, const std::string &name=""){
    dims.resize(1);
    dims[0]=i;  
    recalculate_total_size(name);
  }
  ProductSpace(){
  }
};


class ProductSpaceIndex{
private: //keep private to ensure consistency!
  int I; //big (collective) index
  std::vector<int> v;  //small indices
  int cur_l;  //current subspace whose index is modified by ++ operator
  ProductSpace space;  

public: 
  size_t get_dims_size()const{
    return space.get_dims_size();
  }
 
  int get_I()const{
    return I;
  }

  const std::vector<int> & get_v()const{
    return v;
  }
  const int & operator[] (int l)const {
    check_bounds(l, v.size(),"ProductSpaceIndex::operator[](int l): l");
    return v[l];
  }
 
  void calculate_v_from_I(){
    int J=I;
    for(int l=(int)space.get_dims_size()-1; l>=0; --l){
      v[l]=J%space.get_dim(l);
      J/=space.get_dim(l);
    }
  }
  void calculate_I_from_v(){
    int I=0;
    for(int l=0; l<(int)space.get_dims_size(); ++l){
      I*=space.get_dim(l);
      I+=v[l];
    }
  }

  void remove_space(int l){
    check_bounds(l, space.get_dims_size(), "ProductSpaceIndex::remove_space");
    int sub=1.;
    for(int k=(int)space.get_dims_size()-1; k>l; --k){
      sub*=space.get_dim(k);
    }
    int J= (I/sub/space.get_dim(l));
    J*=sub;
    J+=(I%sub);
    I=J;
    v.erase(v.begin()+l);
    space.erase(l);
  }

  ProductSpaceIndex &operator++(){
    ++I;
//    calculate_v_from_I();

    ++v[cur_l];
    if(v[cur_l]>=space.get_dim(cur_l)){
      v[cur_l]=0;
      if(cur_l<1)return *this; //done (will be signaled by I)
      --cur_l;
      --I;
      return operator++();
    }
    cur_l=(int)space.get_dims_size()-1;

    return *this;
  } 
  bool done()const{
    return (I >= space.get_total_size() );
  }
  void reset(const ProductSpace &space_){
    space=space_;
    space.check_validity("ProductSpace for Index");
    I=0; 
    v.clear();
    v.resize(space.get_dims_size(),0);
    cur_l=(int)v.size()-1;
  }
  void set_from_I(int I_){
    I=I_;
    calculate_v_from_I();
    cur_l=(int)v.size()-1;
  } 
  void set_from_v(std::vector<int> v_){
    v=v_;
    if(v.size()!=space.get_dims_size()){
      std::stringstream ss;
      ss<<"ProductSpaceIndex::set_from_v: v.size()!=dims.size()!"<<std::endl; 
      throw(std::runtime_error(ss.str()));
    }
    calculate_I_from_v();
    cur_l=(int)v.size()-1;
  }
 
  ProductSpaceIndex(const ProductSpace &space){
    reset(space);
  }
};


#endif
