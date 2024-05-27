#ifndef ABSTRACT_TENSOR_DEFINED_H
#define ABSTRACT_TENSOR_DEFINED_H

#include <vector> 
#include <complex>
#include <iostream>
#include <cstdlib>

namespace ACE{

class Tensor_Index;

/// Abstract class describing a tensor, dense or MPS
template <typename T>class Tensor_ScalarType{
public:
 
  virtual int get_rank() const =0;    //how many indices does the tensor have?
  virtual int get_dim(int i) const =0; //dimension for index number i
  virtual int get_total_size() const;
  
  virtual void resize(const std::vector<int> &list) = 0;

  virtual T operator() (const Tensor_Index &ind) const=0;

  void print_dims(std::ostream &os);
  virtual ~Tensor_ScalarType(){}
};


typedef Tensor_ScalarType<std::complex<double> > Tensor;
typedef Tensor_ScalarType<double> Tensor_real;

class Tensor_Index{
public:
  std::vector<int> list;

  inline void resize(size_t s){ list.resize(s);}
  inline int get_rank() const { return list.size(); }
  
  inline int & operator[](size_t i){ return list[i];}
  inline int operator[](size_t i)const{ return list[i];}

  inline int & back(){ return list.back();}
  inline int back()const{ return list.back();}

  bool done(const Tensor &ten);
  void set_done(const Tensor &ten);
  void reset(const Tensor &ten);
  void increment(const Tensor &ten);
  void increment_from_back(const Tensor &ten);
  
  inline Tensor_Index(const std::vector<int> &list_) : list(list_){
  }
  inline Tensor_Index(const Tensor &ten){
    reset(ten);
  }
  inline Tensor_Index(){
  }
  static Tensor_Index from_block(int b, int rank, int dim);
};

std::ostream & operator<<(std::ostream &os, const Tensor_Index & ti);



class Tensor_Dimensions: public Tensor{
protected:
  std::vector<int> dim;
public:

  inline virtual int get_rank() const { return dim.size(); }
  inline virtual int get_dim(int i) const { return dim[i]; }
  inline const std::vector<int> &get_dims() const { return dim; }

  virtual void resize(const std::vector<int> &list);
  
  virtual void resize(int rank_, int dim_);
  
  virtual std::complex<double> operator() (const Tensor_Index &ind)const;

  int get_block_index(const Tensor_Index &ind) const;

  inline Tensor_Dimensions(){}
  inline Tensor_Dimensions(const std::vector<int> & dims_) : dim(dims_) {}
  inline Tensor_Dimensions(int rank_, int dim_){
    resize(rank_, dim_);
  }
  virtual ~Tensor_Dimensions(){}
};
}//namespace

#endif
