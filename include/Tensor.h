#ifndef ABSTRACT_TENSOR_DEFINED_H
#define ABSTRACT_TENSOR_DEFINED_H

#include <vector> 
#include <complex>
#include <iostream>
#include <cstdlib>


class Tensor_Index;

/// Abstract class describing a tensor, dense or MPS
class Tensor{
public:
 
  virtual int get_rank() const =0;    //how many indices does the tensor have?
  virtual int get_dim(int i) const =0; //dimension for index number i
  virtual int get_total_size() const{
    if(get_rank()<1)return 0;
    int b=1; for(size_t i=0; i<get_rank(); i++) b*=get_dim(i);
    return b;
  }
  virtual void resize(const std::vector<int> &list) = 0;

  virtual std::complex<double> operator() (const Tensor_Index &ind) const=0;

  void print_dims(std::ostream &os){
    if(get_rank()<1)return;
    os<<get_dim(0);
    for(size_t i=1; i<get_rank(); i++){
      os<<" "<<get_dim(i);
    }
  }
  virtual ~Tensor(){}
};

class Tensor_Index{
public:
  std::vector<int> list;

  void resize(size_t s){ list.resize(s);}
  int get_rank() const { return list.size(); }
  
  int & operator[](size_t i){ return list[i];}
  int operator[](size_t i)const{ return list[i];}

  int & back(){ return list.back();}
  int back()const{ return list.back();}

  bool done(const Tensor &ten){ return list.back()>=ten.get_dim(ten.get_rank()-1); }
  void set_done(const Tensor &ten){list.back()=ten.get_dim(ten.get_rank()-1);}

  void reset(const Tensor &ten){ list.clear(); list.resize(ten.get_rank(),0); }
  void increment(const Tensor &ten){
    for(size_t i=0; i<list.size(); i++){
      list[i]++; 
      if(list[i]>=ten.get_dim(i) && i<list.size()-1){
        list[i]=0;
      }else return;
    }
  }
  void increment_from_back(const Tensor &ten){
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
  Tensor_Index(const std::vector<int> &list_) : list(list_){
  }
  Tensor_Index(const Tensor &ten){
    reset(ten);
  }
  Tensor_Index(){
  }
  static Tensor_Index from_block(int b, int rank, int dim){
     Tensor_Index index;
     index.resize(rank);
     for(int i=rank-1; i>=0 ; i--){
       index[i]=b%dim;
       b/=dim;
     }
     return index;
  }
};

std::ostream & operator<<(std::ostream &os, const Tensor_Index & ti){
  if(ti.get_rank()<1)return os;
  os<<ti[0];
  for(size_t i=1; i<ti.get_rank(); i++){ os<<" "<<ti[i];}
  return os;
}



class Tensor_Dimensions: public Tensor{
protected:
  std::vector<int> dim;
  std::complex<double> *mem;
public:

  virtual int get_rank() const { return dim.size(); }
  virtual int get_dim(int i) const { return dim[i]; }
  const std::vector<int> &get_dims() const { return dim; }
  virtual void resize(const std::vector<int> &list){
    dim=list;
  }
  virtual void resize(int rank_, int dim_){
    dim.clear();
    dim.resize(rank_,dim_);
  }
  virtual std::complex<double> operator() (const Tensor_Index &ind) const{
    std::cerr<<"ERROR: Trying to access elements of Tensor_Dimensions as Tensor!"<<std::endl;
    exit(1);
  }

  int get_block_index(const Tensor_Index &ind) const{
    if(ind.get_rank()!=dim.size()){
      std::cerr<<"Tensor_Dense: operator(): ind.get_rank()!=dim.size()!"<<std::endl;
      exit(1);
    }
    if(dim.size()<1)return 0.;
    int b=ind[0];
    for(int n=1; n<dim.size(); n++){
      b*=dim[n];
      b+=ind[n];
    }
    return b;
  }

  Tensor_Dimensions(){}
  Tensor_Dimensions(int rank_, int dim_){
    resize(rank_, dim_);
  }
  virtual ~Tensor_Dimensions(){}
};

#include "Tensor_Dense.h"
#include "MPS.h"

#endif
