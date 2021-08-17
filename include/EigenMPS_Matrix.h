#ifndef ACE_EIGEN_MPS_MATRIX_DEFINED_H
#define ACE_EIGEN_MPS_MATRIX_DEFINED_H

#include <Eigen/Dense>
#include <vector>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <fstream>

class EigenMPS_Matrix {
public:
  std::vector<Eigen::MatrixXcd> a;

  int get_dim_o()const{
    return a.size();
  }
  int get_dim_d1()const{
    if(a.size()<1)return 0;
    return a[0].rows();
  }
  int get_dim_d2()const{
    if(a.size()<1)return 0;
    return a[0].cols();
  }
 
  void resize(int o, int d1, int d2){
    if(o<0||d1<0||d2<0){
      std::cerr<<"EigenMPS_Matrix::resize: o<0||d1<0||d2<0: "<<o<<", "<<d1<<", "<<d2<<std::endl;
      exit(1);
    }
    std::vector<Eigen::MatrixXcd> tmp(o, Eigen::MatrixXcd::Zero(d1, d2));
    a.swap(tmp);
  }
  void fill(std::complex<double> arg){
    for(int o=0; o<get_dim_o(); o++){
      for(int d1=0; d1<get_dim_d1(); d1++){
        for(int d2=0; d2<get_dim_d2(); d2++){
          a[o](d1,d2)=arg;
        }
      }
    }
  }

  const Eigen::MatrixXcd & operator[](int o)const{
    return a[o];
  }
  Eigen::MatrixXcd & operator[](int o){
    return a[o];
  }

  std::complex<double> get_check(int o, int d1, int d2){
    if(o<0||o>=a.size()){
      std::cerr<<"EigenMPS_Matrix: o<0||o>=a.size(): "<<o<<"/"<<a.size()<<"!"<<std::endl;
      exit(1);
    }
    if(d1<0||d1>=a[o].rows()){
      std::cerr<<"EigenMPS_Matrix:  d1<0||d1>=a["<<o<<"].rows(): "<<d1<<"/"<<a[o].rows()<<"!"<<std::endl;
      exit(1);
    }
    if(d2<0||d2>=a[o].cols()){
      std::cerr<<"EigenMPS_Matrix:  d2<0||d2>=a["<<o<<"].cols(): "<<d2<<"/"<<a[o].cols()<<"!"<<std::endl;
      exit(1);
    }
  }

  double max_element_abs()const{
    double max=0;
    for(int o=0; o<get_dim_o(); o++){
      for(int d1=0; d1<get_dim_d1(); d1++){
        for(int d2=0; d2<get_dim_d2(); d2++){
          if(abs(a[o](d1,d2))>max)max=abs(a[o](d1,d2));
        }
      }
    }
    return max; 
  }
  void print_dims(std::ostream &os=std::cout)const{
    os<<get_dim_o()<<" "<<get_dim_d1()<<" "<<get_dim_d2();
  }

  void swap(EigenMPS_Matrix &other){
    a.swap(other.a);
  }
  void swap(std::vector<Eigen::MatrixXcd> &other){
    a.swap(other);
  }

  EigenMPS_Matrix(int o, int d1, int d2){
    resize(o, d1, d2);
  }
  EigenMPS_Matrix(){}
  virtual ~EigenMPS_Matrix(){}
};

#endif
