#ifndef MPS_MATRIX_DEFINED_H
#define MPS_MATRIX_DEFINED_H

#include <iostream>
#include <Eigen/Dense>
//#include <fstream>
//#include "CheckMatrix.hpp"

namespace ACE{

template <typename T> class MPS_Matrix_ScalarType{
public:
  int dim_d1, dim_i, dim_d2;
  T *mem;

  inline T &operator()(int i, int d1, int d2){
    return mem[(d1*dim_i+i)*dim_d2+d2];
  }
  inline T &operator()(int i, int d1, int d2)const{
    return mem[(d1*dim_i+i)*dim_d2+d2];
  }
  inline void allocate(){
    mem=new T[dim_d1*dim_d2*dim_i];
  }
  inline void deallocate(){
    delete[] mem;
  }
  inline double norm()const{
    double res=0;
    for(int i=0; i<dim_d1*dim_d2*dim_i; i++){
      res+=std::norm(mem[i]);
    }
    return sqrt(res);
  }

  void resize(int dim_i_, int dim_d1_, int dim_d2_);
  
  void fill(const T &c);
  
  void resize_fill_one(int dim_i_, int dim_d1_, int dim_d2_);
  
  void set_zero();
  
  double max_element_abs()const;
  double norm2_d1(int d1)const;

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> get_Matrix_d1i_d2()const;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> get_Matrix_d1_id2()const;
  void set_from_Matrix_d1i_d2(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &M, int dimi);
  void set_from_Matrix_d1_id2(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &M, int dimi);
 

  inline Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0 , Eigen::OuterStride<> >  get_Matrix_d1_d2(int i)const{
    return Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0 , Eigen::OuterStride<> > (mem+i*dim_d2, dim_d1, dim_d2, Eigen::OuterStride<>(dim_i*dim_d2));
  }
  
  void print_dims(std::ostream &os=std::cout)const;
  
  void print_HR(const std::string &fname, double threshold=0)const;
  
  void copy(const MPS_Matrix_ScalarType<T> &other);

  void swap(MPS_Matrix_ScalarType<T> &other);
  
  void inner_multiply_left(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & M);

  void inner_multiply_right(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & M);

  void contract_front();
  
  void contract_back();
  
  void multiply(const MPS_Matrix_ScalarType<T> &other);

  void write_binary(std::ostream &ofs)const;
  void read_binary(std::istream &ifs);

  inline MPS_Matrix_ScalarType &operator=(const MPS_Matrix_ScalarType<T> &other){
    copy(other);
    return *this; 
  }
  inline MPS_Matrix_ScalarType(const MPS_Matrix_ScalarType<T> &other){
    dim_d1=dim_d2=1;
    dim_i=0;
    allocate();
    copy(other);
  }
  inline MPS_Matrix_ScalarType(int dim_i_, int dim_d1_=1, int dim_d2_=1)
   : dim_d1(dim_d1_), dim_i(dim_i_), dim_d2(dim_d2_){
    allocate();
  }
  inline MPS_Matrix_ScalarType(){
    dim_d1=dim_d2=1;
    dim_i=0;
    allocate();
  }
  ~MPS_Matrix_ScalarType(){
    deallocate();
  }
  static MPS_Matrix_ScalarType<T> get_trivial(int dim_i){
    MPS_Matrix_ScalarType<T> ret(dim_i, 1, 1);
    ret.fill(1.);
    return ret;
  }
};

typedef MPS_Matrix_ScalarType<std::complex<double> > MPS_Matrix;
typedef MPS_Matrix_ScalarType<double> MPS_Matrix_real;

template <typename T>
std::ostream &operator<<(std::ostream &os, const MPS_Matrix_ScalarType<T> &a);

}//namespace

#endif
