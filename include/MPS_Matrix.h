#ifndef MPS_MATRIX_DEFINED_H
#define MPS_MATRIX_DEFINED_H

#include <Eigen/Dense>
#include "CheckMatrix.h"

class MPS_Matrix{
public:
  int dim_d1, dim_i, dim_d2;
  std::complex<double> *mem;

  std::complex<double> &operator()(int i, int d1, int d2){
    return mem[(d1*dim_i+i)*dim_d2+d2];
  }
  const std::complex<double> &operator()(int i, int d1, int d2)const{
    return mem[(d1*dim_i+i)*dim_d2+d2];
  }
  void allocate(){
    mem=new std::complex<double>[dim_d1*dim_d2*dim_i];
  }
  void deallocate(){
    delete[] mem;
  }
  void resize(int dim_i_, int dim_d1_, int dim_d2_){
    dim_i=dim_i_;
    dim_d1=dim_d1_;
    dim_d2=dim_d2_;
    deallocate();
    allocate();
  }
  void fill(const std::complex<double> &c){
    for(int i=0; i<dim_i*dim_d1*dim_d2; i++)mem[i]=c;
  }
  void resize_fill_one(int dim_i_, int dim_d1_, int dim_d2_){
    resize(dim_i_, dim_d1_, dim_d2_);
    fill(1.);
  }
  void set_zero(){
    for(int i=0; i<dim_i*dim_d1*dim_d2; i++)mem[i]=0.;
  }
  double max_element_abs()const{
    double a=0;
    int end=dim_d1*dim_d2*dim_i;
    for(int x=0; x<end; x++){
      if(abs(mem[x])>a)a=abs(mem[x]);
    }
    return a; 
  }
  void print_dims(std::ostream &os=std::cout)const{
    os<<dim_i<<" "<<dim_d1<<" "<<dim_d2;
  }
  void print_HR(const std::string &fname)const{
    std::ofstream ofs(fname.c_str());
    ofs<<"Dimensions: "<<dim_i<<" "<<dim_d1<<" "<<dim_d2<<std::endl;
    for(int i=0; i<dim_i; i++){
      for(int d1=0; d1<dim_d1; d1++){
        for(int d2=0; d2<dim_d2; d2++){
          ofs<<dim_i<<" "<<dim_d1<<" "<<dim_d2<<" "<<operator()(i,d1,d2)<<std::endl;
        }
      }
    }
  }
  void copy(const MPS_Matrix &other){
    resize(other.dim_i, other.dim_d1, other.dim_d2);
    for(int i=0; i<dim_i*dim_d1*dim_d2; i++){
      mem[i]=other.mem[i];
    }
  }

  void swap(MPS_Matrix &other){
    int dim_i_=dim_i;
    int dim_d1_=dim_d1;
    int dim_d2_=dim_d2;
    std::complex<double> *mem_=mem;
    dim_i=other.dim_i;
    dim_d1=other.dim_d1;
    dim_d2=other.dim_d2;
    mem=other.mem;
    other.dim_i=dim_i_;
    other.dim_d1=dim_d1_;
    other.dim_d2=dim_d2_;
    other.mem=mem_;
  }
  void inner_multiply_left(const Eigen::MatrixXcd & M){
    if(dim_d1!=M.cols()){
      std::cerr<<"MPS_Matrix::inner_multiply_left: dim_d1!=M.cols()!"<<std::endl;
      exit(1);
    }
    MPS_Matrix A(dim_i, M.rows(), dim_d2); 
    for(int i=0; i<dim_i; i++){
      for(int r=0; r<M.rows(); r++){
        for(int d2=0; d2<dim_d2; d2++){
          A(i,r,d2)=0.;
          for(int d1=0; d1<dim_d1; d1++){
            A(i,r,d2)+=M(r,d1)*operator()(i,d1,d2);
          }
        }
      }
    }
    swap(A);
  }
  void inner_multiply_right(const Eigen::MatrixXcd & M){
    if(dim_d2!=M.rows()){
      std::cerr<<"MPS_Matrix::inner_multiply_right: dim_d2!=M.rows()!"<<std::endl;
      exit(1);
    }
    MPS_Matrix A(dim_i, dim_d1, M.cols()); 
    for(int i=0; i<dim_i; i++){
      for(int d1=0; d1<dim_d1; d1++){
        for(int c=0; c<M.cols(); c++){
          A(i,d1,c)=0.;
          for(int d2=0; d2<dim_d2; d2++){
            A(i,d1,c)+=operator()(i,d1,d2)*M(d2,c);
          }
        }
      }
    }
    swap(A);
  }
  void contract_front(){
    if(dim_d1<2)return;
    MPS_Matrix A(dim_i, 1, dim_d2);
    A.set_zero();
    for(int i=0; i<dim_i; i++){
      for(int d1=0; d1<dim_d1; d1++){
        for(int d2=0; d2<dim_d2; d2++){
          A(i,0,d2)+=operator()(i,d1,d2);
        }
      }
    }
    swap(A);
  }
  void contract_back(){
    if(dim_d2<2)return;
    MPS_Matrix A(dim_i, dim_d1, 1);
    A.set_zero();
    for(int i=0; i<dim_i; i++){
      for(int d1=0; d1<dim_d1; d1++){
        for(int d2=0; d2<dim_d2; d2++){
          A(i,d1,0)+=operator()(i,d1,d2);
        }
      }
    }
    swap(A);
  }
  void multiply(const MPS_Matrix &other){
    if(other.dim_i!=dim_i){
      std::cerr<<"MPS_Matrix::multiply: other.dim_i!=dim_i!"<<std::endl;
      exit(1);
    }

    MPS_Matrix A(dim_i, dim_d1*other.dim_d1, dim_d2*other.dim_d2);
#ifdef PRINT_MPS_MATRIX_MULTIPLY_DIMS
std::cout<<"MPS_Matrix::multiply: dims: "<<A.dim_i<<" "<<A.dim_d1<<" "<<A.dim_d2<<" <- "<<dim_d1<<"*"<<other.dim_d1<<", "<<dim_d2<<"*"<<other.dim_d2<<std::endl;
#endif
    //A.set_zero();
    for(int i=0; i<dim_i; i++){
      for(int d1=0; d1<dim_d1; d1++){
        for(int d2=0; d2<dim_d2; d2++){
          for(int od1=0; od1<other.dim_d1; od1++){
            for(int od2=0; od2<other.dim_d2; od2++){
              A(i, d1*other.dim_d1 + od1, d2*other.dim_d2 + od2)=
                 operator()(i,d1,d2) * other(i,od1,od2);
            }
          }
        }
      }
    }
    swap(A);
  }

  MPS_Matrix &operator=(const MPS_Matrix &other){
    copy(other);
    return *this; 
  }
  MPS_Matrix(const MPS_Matrix &other){
    dim_d1=dim_d2=1;
    dim_i=0;
    allocate();
    copy(other);
  }
  MPS_Matrix(int dim_i_, int dim_d1_=1, int dim_d2_=1)
   : dim_d1(dim_d1_), dim_i(dim_i_), dim_d2(dim_d2_){
    allocate();
  }
  MPS_Matrix(){
    dim_d1=dim_d2=1;
    dim_i=0;
    allocate();
  }
  ~MPS_Matrix(){
    deallocate();
  }
};
std::ostream &operator<<(std::ostream &os, const MPS_Matrix &a){
  a.print_dims(os);
  return os;
}



#endif
