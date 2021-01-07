#ifndef TENSOR_DENSE_DEFINED_H
#define TENSOR_DENSE_DEFINED_H

#include "Tensor.h"

class Tensor_Dense: public Tensor{
protected:
  std::vector<int> dim;
  std::complex<double> *mem;

public:
  virtual int get_rank() const { return dim.size(); }
  virtual int get_dim(int i) const { return dim[i]; }
  const std::vector<int> &get_dims() const { return dim; }

  /* backwards:
  int get_block_index(const Tensor_Index &ind) const{
    if(ind.get_rank()!=dim.size()){
      std::cerr<<"Tensor_Dense: operator(): ind.get_rank()!=dim.size()!"<<std::endl;
      exit(1);
    } 
    if(dim.size()<1)return 0.;
    int b=ind[dim.size()-1];
    for(int n=(int)dim.size()-2; n>=0; n--){
      b*=dim[n];
      b+=ind[n];
    }
    return b;
  } */

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

  virtual std::complex<double> operator() (const Tensor_Index &ind) const{
    int b=get_block_index(ind);
    return mem[b];
  }
  //In contrast to, e.g., MPS versions: write access to elements
  virtual std::complex<double> &operator[] (int b){
    return mem[b];
  }
  virtual const std::complex<double> &operator[] (int b) const{
    return mem[b];
  }
  std::complex<double> &operator[](const Tensor_Index &ind){
    int b=get_block_index(ind);
    return mem[b];
  }

  void allocate(){ mem=new std::complex<double>[get_total_size()]; }
  void deallocate(){ if(mem!=NULL)delete[] mem; mem=NULL;}
  virtual void resize(const std::vector<int> &list){
    deallocate();
    dim=list;
    allocate();
  }
  virtual void resize(int rank_, int dim_){
    deallocate();
    dim.clear();
    dim.resize(rank_,dim_);
    allocate();
  }
  void fill(std::complex<double> c){
    for(int i=0; i<get_total_size(); i++)mem[i]=c;
  }
  void copy(const Tensor_Dense &other){
    resize(other.dim);
    for(int i=0; i<get_total_size(); i++)mem[i]=other.mem[i];
  }
  Tensor_Dense &operator=(const Tensor_Dense &other){
    copy(other);
    return *this;
  }

  void swap(Tensor_Dense &other){
    if(other.dim.size()!=dim.size()){
      std::cerr<<"Tensor_Dense::swap: other.dim.size()!=dim.size()!"<<std::endl;
      exit(1);
    }
    for(size_t i=0; i<dim.size(); i++){
      if(other.dim[i]!=dim[i]){
        std::cerr<<"Tensor_Dense::swap: other.dim[i]!=dim[i] for i="<<i<<std::endl;
        exit(1);
      }
    }
    std::complex<double> *swp=mem;
    mem=other.mem; 
    other.mem=swp;
  }

  Tensor_Dense(const std::vector<int> &v) :dim(v){
    allocate();
  }
  Tensor_Dense(int rank, int dim){
    allocate();
    resize(rank,dim);
  }
  Tensor_Dense(const Tensor_Dense &other){
    allocate();
    copy(other);
  }
  Tensor_Dense(){
    allocate();
  }
  virtual ~Tensor_Dense(){
    deallocate();
  }

};



#endif
