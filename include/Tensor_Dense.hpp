#ifndef TENSOR_DENSE_DEFINED_H
#define TENSOR_DENSE_DEFINED_H

#include "Tensor.hpp"
//#include <vector>
//#include <complex>

namespace ACE{

class Tensor_Dense: public Tensor{
protected:
  std::vector<int> dim;
  std::complex<double> *mem;

public:
  inline virtual int get_rank() const { return dim.size(); }
  inline virtual int get_dim(int i) const { return dim[i]; }
  inline const std::vector<int> &get_dims() const { return dim; }

  // backwards:
  //int get_block_index(const Tensor_Index &ind)const; 

  int get_block_index(const Tensor_Index &ind)const;

  virtual std::complex<double> operator() (const Tensor_Index &ind)const;
  
  //In contrast to, e.g., MPS versions: write access to elements
  inline virtual std::complex<double> &operator[] (int b){
    return mem[b];
  }
  virtual const std::complex<double> &operator[] (int b) const{
    return mem[b];
  }
  std::complex<double> &operator[](const Tensor_Index &ind);

  void allocate();
  void deallocate();
  virtual void resize(const std::vector<int> &list);
  virtual void resize(int rank_, int dim_);
  void fill(std::complex<double> c);
  void copy(const Tensor_Dense &other);

  inline Tensor_Dense &operator=(const Tensor_Dense &other){
    copy(other);
    return *this;
  }

  void swap(Tensor_Dense &other);

  inline Tensor_Dense(const std::vector<int> &v) :dim(v){
    allocate();
  }
  inline Tensor_Dense(int rank, int dim){
    allocate();
    resize(rank,dim);
  }
  inline Tensor_Dense(const Tensor_Dense &other){
    allocate();
    copy(other);
  }
  inline Tensor_Dense(){
    allocate();
  }
  inline virtual ~Tensor_Dense(){
    deallocate();
  }

};


}//namespace
#endif
