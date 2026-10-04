#ifndef FOURLEG_DEFINED_H
#define FOURLEG_DEFINED_H

#include <iostream>
#include <Eigen/Dense>
#include "MPS_Matrix.hpp"

namespace ACE{

template <typename T> class FourLeg_ScalarType{
public:
  //   dim3 |
  //dim2 --   -- dim4
  //   dim1 |
  int dim_d1, dim_d2, dim_d3, dim_d4;
  T *mem;

  inline T &operator()(int d1, int d2, int d3, int d4){
    return mem[((d1*dim_d2+d2)*dim_d3+d3)*dim_d4+d4];
  }
  inline T &operator()(int d1, int d2, int d3, int d4)const{
    return mem[((d1*dim_d2+d2)*dim_d3+d3)*dim_d4+d4];
  }
  inline void allocate(){
    mem=new T[dim_d1*dim_d2*dim_d3*dim_d4];
  }
  inline void deallocate(){
    delete[] mem;
  }
  inline double norm()const{
    double res=0;
    for(int i=0; i<dim_d1*dim_d2*dim_d3*dim_d4; i++){
      res+=std::norm(mem[i]);
    }
    return sqrt(res);
  }
  inline void scale(T c){
    for(int i=0; i<dim_d1*dim_d2*dim_d3*dim_d4; i++){
      mem[i]*=c;
    }
  }

  void resize(int dim_d1_, int dim_d2_, int dim_d3_, int dim_d4_);
  
  void fill(const T &c);
  
  void resize_fill_one(int dim_d1_, int dim_d2_, int dim_d3_, int dim_d4_);
  
  void set_zero();
  
  double max_element_abs()const;

//  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> get_Matrix_d1i_d2()const;
//  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> get_Matrix_d1_id2()const;
//  void set_from_Matrix_d1i_d2(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &M, int dimi);
//  void set_from_Matrix_d1_id2(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &M, int dimi);
 

//  inline Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0 , Eigen::OuterStride<> >  get_Matrix_d1_d2(int i)const{
//    return Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0 , Eigen::OuterStride<> > (mem+i*dim_d2, dim_d1, dim_d2, Eigen::OuterStride<>(dim_i*dim_d2));
//  }
  
  void print_dims(std::ostream &os=std::cout)const;
  
  void print_HR(const std::string &fname, double threshold=0)const;
  
  void copy(const FourLeg_ScalarType<T> &other);
  //resize and copy truncated version of "other" or zero-pad.
  void resize_copy(int new_d1, int new_d2, int new_d3, int new_d4, const FourLeg_ScalarType<T> &other);

  void swap(FourLeg_ScalarType<T> &other);
 
  //assume: contraction over first index
  void multiply_d1(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & M);
  void multiply_d2(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & M);
  void multiply_d3(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & M);
  void multiply_d4(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & M);

  void join_d3_d1(const FourLeg_ScalarType<T> &other);
  void join_d3_d1_multiply_d2(const FourLeg_ScalarType<T> &other, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & P);

  //Reduce to three-leg object by fixing d1; mapping d3 -> i; d2 -> d1; d4 -> d3;
  MPS_Matrix_ScalarType<T> reduce_fix_d1(int d1) const;

  //cast to Matrix  
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> get_Matrix_d1_d2d3d4()const;
  void set_from_Matrix_d1_d2d3d4(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &M, int d2, int d3);

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> get_Matrix_d1d2d3_d4()const;
  void set_from_Matrix_d1d2d3_d4(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &M, int d2, int d3);

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> get_Matrix_d2_d1d3d4()const;
  void set_from_Matrix_d2_d1d3d4(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &M, int d3, int d4);

  Eigen::VectorXcd get_reduced_d3(const Eigen::VectorXcd &cap_d1, const Eigen::VectorXcd &cap_d2, const Eigen::VectorXcd &cap_d4)const;


  void write_binary(std::ostream &ofs)const;
  void read_binary(std::istream &ifs);

  inline FourLeg_ScalarType &operator=(const FourLeg_ScalarType<T> &other){
    copy(other);
    return *this; 
  }
  inline FourLeg_ScalarType(const FourLeg_ScalarType<T> &other){
    dim_d1=0; dim_d2=dim_d3=dim_d4=1;
    allocate();
    copy(other);
  }
  inline FourLeg_ScalarType(int dim_d1_, int dim_d2_=1, int dim_d3_=1, int dim_d4_=1)
   : dim_d1(dim_d1_), dim_d2(dim_d2_), dim_d3(dim_d3_), dim_d4(dim_d4_){
    allocate();
  }
  inline FourLeg_ScalarType(){
    dim_d1=0; dim_d2=dim_d3=dim_d4=1;
    allocate();
  }
  ~FourLeg_ScalarType(){
    deallocate();
  }
  static FourLeg_ScalarType<T> get_trivial(){
    FourLeg_ScalarType<T> ret(1, 1, 1, 1);
    ret.fill(1.);
    return ret;
  }
  static FourLeg_ScalarType<T> Zero(int d1, int d2, int d3, int d4){
    FourLeg_ScalarType<T> ret(d1, d2, d3, d4);
    ret.set_zero();
    return ret;
  }
  static FourLeg_ScalarType<T> Identity_d1_d3(int d);

  static FourLeg_ScalarType<T> Operator_d1_d3(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & M);
  static FourLeg_ScalarType<T> Operator_forward_d1_d3(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & M);
  static FourLeg_ScalarType<T> Operator_backward_d1_d3(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & M);
  static FourLeg_ScalarType<T> Operator_forward_backward_d1_d3(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & op1, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & op2);
};

typedef FourLeg_ScalarType<std::complex<double> > FourLeg;
typedef FourLeg_ScalarType<double> FourLeg_real;

template <typename T> std::ostream &operator<<(std::ostream &os, const FourLeg_ScalarType<T> &a);


FourLeg MPS_Matrix_to_FourLeg(const MPS_Matrix &M, int dim_d2);
MPS_Matrix FourLeg_to_MPS_Matrix(const FourLeg &f);
}//namespace

#endif
