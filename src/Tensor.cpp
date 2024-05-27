#include "Tensor.hpp"
#include <vector> 
#include <complex>
#include <iostream>
#include <cstdlib>

namespace ACE{

/// Abstract class describing a tensor, dense or MPS
template <typename T>
  int Tensor_ScalarType<T>::get_total_size() const{
    if(get_rank()<1)return 0;
    int b=1; for(int i=0; i<get_rank(); i++) b*=get_dim(i);
    return b;
  }

template <typename T>
  void Tensor_ScalarType<T>::print_dims(std::ostream &os){
    if(get_rank()<1)return;
    os<<get_dim(0);
    for(int i=1; i<get_rank(); i++){
      os<<" "<<get_dim(i);
    }
  }


//Tensor_Index:


  bool Tensor_Index::done(const Tensor &ten){ 
    return list.back()>=ten.get_dim(ten.get_rank()-1); 
  }
  void Tensor_Index::set_done(const Tensor &ten){
    list.back()=ten.get_dim(ten.get_rank()-1);
  }
  void Tensor_Index::reset(const Tensor &ten){ 
    list.clear(); 
    list.resize(ten.get_rank(),0); 
  }

  void Tensor_Index::increment(const Tensor &ten){
    for(size_t i=0; i<list.size(); i++){
      list[i]++; 
      if(list[i]>=ten.get_dim(i) && i<list.size()-1){
        list[i]=0;
      }else return;
    }
  }

  void Tensor_Index::increment_from_back(const Tensor &ten){
    if(done(ten))return;
    for(int i=(int)list.size()-1; i>=0; i--){
      list[i]++; 
      if(list[i]>=ten.get_dim(i)){
        if(i==0){
          set_done(ten);
          return;
        }
        list[i]=0;
      }else return;
    }
  }
  Tensor_Index Tensor_Index::from_block(int b, int rank, int dim){
     Tensor_Index index;
     index.resize(rank);
     for(int i=rank-1; i>=0 ; i--){
       index[i]=b%dim;
       b/=dim;
     }
     return index;
  }

std::ostream & operator<<(std::ostream &os, const Tensor_Index & ti){
  if(ti.get_rank()<1)return os;
  os<<ti[0];
  for(int i=1; i<ti.get_rank(); i++){ os<<" "<<ti[i];}
  return os;
}


// Tensor_Dimensions: 

  void Tensor_Dimensions::resize(const std::vector<int> &list){
    dim=list;
  }
  void Tensor_Dimensions::resize(int rank_, int dim_){
    dim.clear();
    dim.resize(rank_,dim_);
  }
  std::complex<double> Tensor_Dimensions::operator() (const Tensor_Index &ind) const{
    std::cerr<<"ERROR: Trying to access elements of Tensor_Dimensions as Tensor!"<<std::endl;
    exit(1);
  }

  int Tensor_Dimensions::get_block_index(const Tensor_Index &ind) const{
    if(ind.get_rank()!=(int)dim.size()){
      std::cerr<<"Tensor_Dense: operator(): ind.get_rank()!=dim.size()!"<<std::endl;
      exit(1);
    }
    if(dim.size()<1)return 0.;
    int b=ind[0];
    for(size_t n=1; n<dim.size(); n++){
      b*=dim[n];
      b+=ind[n];
    }
    return b;
  }

template class Tensor_ScalarType<double>;
template class Tensor_ScalarType<std::complex<double>>;
}//namespace
