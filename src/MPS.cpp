#include "MPS.hpp"
#include "Tensor.hpp"
#include <Eigen/Core>
#include <Eigen/SVD>
#include <iostream>
#include <fstream>
#include <complex>
#include "MPS_Matrix.hpp"
#include "TimeGrid.hpp"
#include "CheckMatrix.hpp"
#include "Compress_Trafo_At.hpp"
#include "Sweep_Trafo_Processor.hpp"

namespace ACE{

template <typename T> 
  void MPS_ScalarType<T>::resize(const std::vector<int> &list){
    a.clear();
    a.resize(list.size());
    for(size_t i=0; i<list.size(); i++){
      a[i].resize(list[i],1,1);
    }
  }
template <typename T> 
  void MPS_ScalarType<T>::resize(int rank, int dim){
    a.clear();
    a.resize(rank);
    for(int i=0; i<rank; i++){
      a[i].resize(dim,1,1);
    }
  }
template <typename T> 
  void MPS_ScalarType<T>::resize_fill_one(int rank, int dim){
    a.clear();
    a.resize(rank);
    for(int i=0; i<rank; i++){
      a[i].resize(dim,1,1);
      a[i].fill(1.);
    }
  }
template <typename T> 
  void MPS_ScalarType<T>::reduce_last(const Eigen::Matrix<T, Eigen::Dynamic, 1> &v){
    if(a.size()<1)return;
    MPS_Matrix_ScalarType<T> &a_ref=a.back();

    MPS_Matrix_ScalarType<T> A(a_ref.dim_i, a_ref.dim_d1, 1);
    check_vector_dim_eq(v, a_ref.dim_d2, "MPS::reduce_last: v");
    A.set_zero();
    for(int i=0; i<a_ref.dim_i; i++){
      for(int d1=0; d1<a_ref.dim_d1; d1++){
        for(int d2=0; d2<a_ref.dim_d2; d2++){
          A(i, d1, 0)+=a_ref(i, d1, d2)*v(d2);
        }
      }
    }
    a_ref.swap(A);
  }
template <typename T> 
  void MPS_ScalarType<T>::reduce_first(const Eigen::Matrix<T, Eigen::Dynamic, 1> &v){
    if(a.size()<1)return;
    MPS_Matrix_ScalarType<T> &a_ref=a[0];

    MPS_Matrix_ScalarType<T> A(a_ref.dim_i, 1, a_ref.dim_d2);
    check_vector_dim_eq(v, a_ref.dim_d1, "MPS::reduce_first: v");
    A.set_zero();
    for(int i=0; i<a_ref.dim_i; i++){
      for(int d1=0; d1<a_ref.dim_d1; d1++){
        for(int d2=0; d2<a_ref.dim_d2; d2++){
          A(i, 0, d2)+=a_ref(i, d1, d2)*v(d1);
        }
      }
    }
    a_ref.swap(A);
  }

template <typename T> 
  void MPS_ScalarType<T>::print_dims(std::ostream &os) const{
    for(size_t i=0; i<a.size(); i++){
      a[i].print_dims(os);
      os<<std::endl;
    }
  }
template <typename T> 
  void MPS_ScalarType<T>::print_dims(const std::string &filename) const{
    std::ofstream ofs(filename.c_str());
    print_dims(ofs);
  }
template <typename T> 
  int MPS_ScalarType<T>::get_max_dim()const{
    int max_dim=0;
    for(size_t i=0; i<a.size(); i++){
      if(a[i].dim_d2>max_dim)max_dim=a[i].dim_d2;
    }
    return max_dim;
  }
template <typename T> 
  int MPS_ScalarType<T>::get_max_dim_pos()const{
    int max_dim=0;
    int pos=-1;
    for(size_t i=0; i<a.size(); i++){
      if(a[i].dim_d2>max_dim){
        max_dim=a[i].dim_d2;
        pos=i;
      }
    }
    return pos;
  }
template <typename T> 
  void MPS_ScalarType<T>::print_max_dim(std::ostream &os) const{
    os<<get_max_dim();
  }
template <typename T> 
  void MPS_ScalarType<T>::print_inner_dims(std::ostream &os) const{
    for(size_t i=0; i<a.size(); i++){
      if(i==0){
        os<<a[i].dim_d1<<" "<<a[i].dim_d2;
      }else{
        os<<" "<<a[i].dim_d2;
      }
    }
  }
template <typename T> 
  double MPS_ScalarType<T>::max_element_abs()const{
    double max=0;
    for(size_t i=0; i<a.size(); i++){
      double m=a[i].max_element_abs();
      if(m>max)max=m;
    }
    return max;
  }

template <typename T> 
  void MPS_ScalarType<T>::check_consistency(const std::vector<int> &i_list)const{
    if(a.size()>0)if(a[0].dim_d1!=1){
      std::cerr<<"MPS::check_consistency: a[0].dim_d1!=1"<<std::endl;
      exit(1);
    }
    if(a.size()>0)if(a[a.size()-1].dim_d2!=1){
      std::cerr<<"MPS::check_consistency: a[last].dim_d2!=1"<<std::endl;
      exit(1);
    }
    for(int k=0; k<((int)a.size())-1; k++){
      if(a[k].dim_d2!=a[k+1].dim_d1){
        std::cerr<<"MPS::check_consistency: a[k].dim_d2!=a[k+1].dim_d1 for k="<<k<<"!"<<std::endl;
        exit(1);
      }
    }
    if(i_list.size()>0){
      if(i_list.size()!=a.size()){
        std::cerr<<"MPS::check_consistency: i_list.size()!=a.size() !"<<std::endl;
        exit(1);
      }
      for(size_t k=0; k<a.size(); k++){
        if(i_list[k]>=a[k].dim_i||i_list[k]<0){
          std::cerr<<"MPS::check_consistency: i_list: out of bounds!"<<std::endl;
          exit(1);
        }
      }
    }
  }

template <typename T> 
  T MPS_ScalarType<T>::operator()(const std::vector<int> &i_list)const{
    check_consistency(i_list);
    std::vector<T> v(1,1.);
    for(size_t k=0; k<a.size(); k++){
      std::vector<T> w(a[k].dim_d2, 0.);
      for(int d2=0; d2<a[k].dim_d2; d2++){
        for(int d1=0; d1<a[k].dim_d1; d1++){
          w[d2]+=v[d1]*a[k](i_list[k], d1, d2);
        } 
      }
      w.swap(v);
    }
    return v[0];
  }

template <typename T> 
  void MPS_ScalarType<T>::multiply_front(const MPS_ScalarType<T> & other){
    if(a.size()<other.a.size()){
      std::cerr<<"MPS::multiply_front: a.size()<other.a.size()!"<<std::endl;
      exit(1);
    }
    for(size_t i=0; i<other.a.size(); i++){
      a[i].multiply(other.a[i]);
    }
  }

template <typename T> 
  void MPS_ScalarType<T>::low_end_multiply(const MPS_ScalarType<T> & other){
    size_t sz;
    if(a.size()<other.a.size()){
      sz=a.size();
    }else{
      sz=other.a.size();
    }
    for(size_t i=0; i<sz; i++){
      a[i].multiply(other.a[i]);
    }
  }

template <typename T> 
  void MPS_ScalarType<T>::insert(int pos, const MPS_Matrix_ScalarType<T> &b, int n_insert){
    if(pos<0){
      std::cerr<<"MPS::insert called with pos<0!"<<std::endl;
      exit(1);
    }
    if(n_insert<0){
      std::cerr<<"MPS::insert called with n_insert<0!"<<std::endl;
      exit(1);
    }
    size_t oldsize=a.size();
    std::vector<MPS_Matrix_ScalarType<T> > mvec(oldsize+n_insert);
    for(size_t n=0; n<(size_t)pos; n++){
      mvec[n].swap(a[n]);
    }
    for(size_t n=pos; n<(size_t)(pos+n_insert); n++){
      mvec[n].copy(b);
    }
    for(size_t n=pos; n<oldsize; n++){
      mvec[n_insert+n].swap(a[n]);
    }
    
    a.clear();
    a.resize(oldsize+n_insert);
    for(size_t n=0; n<(size_t)(oldsize+n_insert); n++){
      a[n].swap(mvec[n]);
    }
  }

template <typename T> 
  void MPS_ScalarType<T>::write_binary(std::ostream &ofs)const{
    size_t sz=a.size();
//std::cout<<"sz: "<<sz<<std::endl;
    ofs.write((char*)&sz, sizeof(int));
    for(size_t i=0; i<sz; i++){
      ofs.write((char*)&a[i].dim_i, sizeof(int));
      ofs.write((char*)&a[i].dim_d1, sizeof(int));
      ofs.write((char*)&a[i].dim_d2, sizeof(int));
    }
    for(size_t i=0; i<sz; i++){
      ofs.write((char*)a[i].mem, sizeof(T)*a[i].dim_i*a[i].dim_d1*a[i].dim_d2);
    }
  }
template <typename T> 
  void MPS_ScalarType<T>::write_binary(const std::string &filename)const{
    std::ofstream ofs(filename.c_str(), std::ios::binary);
    write_binary(ofs);
  }

template <typename T> 
  void MPS_ScalarType<T>::read_binary(std::istream &ifs){
    int itmp;
    ifs.read((char*)&itmp, sizeof(int));
    std::cout<<"Reading size: "<<itmp<<std::endl;
    a.resize(itmp);
    for(size_t i=0; i<a.size(); i++){
      int dim_i, dim_d1, dim_d2;
      ifs.read((char*)&dim_i, sizeof(int));
      ifs.read((char*)&dim_d1, sizeof(int));
      ifs.read((char*)&dim_d2, sizeof(int));
      a[i].resize(dim_i, dim_d1, dim_d2);
    }
    for(size_t i=0; i<a.size(); i++){
      ifs.read((char*)a[i].mem, sizeof(T)*a[i].dim_i*a[i].dim_d1*a[i].dim_d2);
    }
  }
template <typename T> 
  void MPS_ScalarType<T>::read_binary(const std::string &filename){
    std::ifstream ifs(filename.c_str(), std::ios::binary);
    read_binary(ifs);
  }
template <typename T> 
  void MPS_ScalarType<T>::copy(const MPS_ScalarType<T> &other){
    a.clear();
    a.resize(other.a.size());
    for(size_t i=0; i<a.size(); i++){
      a[i].copy(other.a[i]);
    }
  }
template <typename T>
  MPS_ScalarType<T> MPS_ScalarType<T>::get_trivial(int nmax, int dim_i){
    MPS_ScalarType<T> ret; 
    ret.a.resize(nmax, MPS_Matrix_ScalarType<T>::get_trivial(dim_i));
    return ret;
  }


template class MPS_ScalarType<double>;
template class MPS_ScalarType<std::complex<double>>;
}//namespace
