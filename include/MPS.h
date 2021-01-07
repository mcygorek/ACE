#ifndef MPS_DEFINED_H_
#define MPS_DEFINED_H_

#include "Tensor.h"
#include <Eigen/Core>
#include <Eigen/SVD>
#include <iostream>
#include <complex>
#include "RankCompressor.h"
#include "MPS_Matrix.h"
#include "IF_TimeGrid.h"
#include "Compress_Trafo_At.h"

/** Goal:

We want to be able to store (and propagate) a tensor
in the form of a matrix-product state:

A_{i_0, i_1, ... i_{n-1}}   = 
\sum_{d_1,d_2,.. d_{n-1}} a_{i_0}^{0 d_1} a_{i_1}^{d_1d_2} .. a_{i_{n-1}}^{d_{n-1} 0 } 

The decomposition will be performed using sweeps of SVDs.


Every "a" is an object with three different indices

*/

class MPS: public Tensor{
public:
  std::vector<MPS_Matrix> a;
 
  virtual int get_rank() const{return a.size();}
  virtual int get_dim(int i) const{return a[i].dim_i;}
  virtual void resize(const std::vector<int> &list){
    a.clear();
    a.resize(list.size());
    for(size_t i=0; i<list.size(); i++){
      a[i].resize(list[i],1,1);
    }
  }
  void resize(int rank, int dim){
    a.clear();
    a.resize(rank);
    for(size_t i=0; i<rank; i++){
      a[i].resize(dim,1,1);
    }
  }
  void resize_fill_one(int rank, int dim){
    a.clear();
    a.resize(rank);
    for(int i=0; i<rank; i++){
      a[i].resize(dim,1,1);
      a[i].fill(1.);
    }
  }

  void print_dims(std::ostream &os=std::cout) const{
    for(size_t i=0; i<a.size(); i++){
      a[i].print_dims(os);
      os<<std::endl;
    }
  }
  int get_max_dim()const{
    int max_dim=0;
    for(size_t i=0; i<a.size(); i++){
      if(a[i].dim_d2>max_dim)max_dim=a[i].dim_d2;
    }
    return max_dim;
  }
  void print_max_dim(std::ostream &os=std::cout) const{
    os<<get_max_dim();
  }

  void check_consistency(const std::vector<int> &i_list=std::vector<int>())const{
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
  std::complex<double> operator()(const std::vector<int> &i_list)const{
    check_consistency(i_list);
    std::vector<std::complex<double> > v(1,1.);
    for(size_t k=0; k<a.size(); k++){
      std::vector<std::complex<double> > w(a[k].dim_d2, 0.);
      for(size_t d2=0; d2<a[k].dim_d2; d2++){
        for(size_t d1=0; d1<a[k].dim_d1; d1++){
          w[d2]+=v[d1]*a[k](i_list[k], d1, d2);
        } 
      }
      w.swap(v);
    }
    return v[0];
  }
  std::complex<double> operator()(const Tensor_Index &ind)const{
    return operator()(ind.list);
  }

  void multiply_front(const MPS & other){
    if(a.size()<other.a.size()){
      std::cerr<<"MPS::multiply_front: a.size()<other.a.size()!"<<std::endl;
      exit(1);
    }
    for(size_t i=0; i<other.a.size(); i++){
      a[i].multiply(other.a[i]);
    }
  }


  void sweep_block_low_to_high(int n, RankCompressor &compressor, double keep_weight=0., Eigen::MatrixXcd *getR=NULL){
#ifdef PRINT_SVD
std::cout<<"sweep low->high: "<<n<<std::endl;
#endif

    if(n<0)return;
    if(n>=a.size()-1)return;

    Eigen::MatrixXcd L, R;
    compressor.compress(a[n], L, R, true, keep_weight);
    int newdim=L.cols();

    //modify a[n]
    {
    MPS_Matrix ar(a[n].dim_i, a[n].dim_d1, newdim);
    for(int d1=0; d1<a[n].dim_d1; d1++){
      for(int i=0; i<a[n].dim_i; i++){
        for(int dn=0; dn<newdim; dn++){
          ar(i,d1,dn) = L(d1*a[n].dim_i+i, dn);
        }
      }
    }
    a[n].swap(ar);   
    }
   
    //modify a[n+1]
    {
    MPS_Matrix ar(a[n+1].dim_i, newdim, a[n+1].dim_d2);
    ar.set_zero();
    for(int d1=0; d1<a[n+1].dim_d1; d1++){
      for(int d2=0; d2<a[n+1].dim_d2; d2++){
        for(int i=0; i<a[n+1].dim_i; i++){
          for(int dn=0; dn<newdim; dn++){
            ar(i,dn,d2) += R(dn, d1) * a[n+1](i, d1, d2);
          }
        }
      }
    }
    a[n+1].swap(ar);   
    }
   
    //apply to environment operators
    if(getR!=NULL){ *getR=R; }
//    if(cta!=NULL){ cta->set_r_if_correct_n(r,n); }
  }


  void sweep_block_high_to_low(int n, RankCompressor &compressor, double keep_weight=0., Eigen::MatrixXcd *getL=NULL){
    if(n<1)return;
#ifdef PRINT_SVD
std::cout<<"sweep high->low: "<<n<<std::endl;
#endif

    Eigen::MatrixXcd L, R;

    compressor.compress(a[n], L, R, false, keep_weight);

    int newdim=L.cols();

    //modify a[n]
    {
    MPS_Matrix ar(a[n].dim_i, newdim, a[n].dim_d2);
    for(int d2=0; d2<a[n].dim_d2; d2++){
      for(int i=0; i<a[n].dim_i; i++){
        for(int dn=0; dn<newdim; dn++){
          ar(i,dn,d2) = R(dn, i*a[n].dim_d2+d2);
        }
      }
    }
    a[n].swap(ar);   
    }
   
    //modify a[n-1]
    {
    MPS_Matrix ar(a[n-1].dim_i, a[n-1].dim_d1, newdim);
    ar.set_zero();
    for(int d1=0; d1<a[n-1].dim_d1; d1++){
      for(int d2=0; d2<a[n-1].dim_d2; d2++){
        for(int i=0; i<a[n-1].dim_i; i++){
          for(int dn=0; dn<newdim; dn++){
            ar(i,d1,dn) += a[n-1](i, d1, d2) * L(d2, dn);
          }
        }
      }
    }
    a[n-1].swap(ar);   
    }

    //apply to environment operators
    if(getL!=NULL){ *getL=L; }
//    if(cta!=NULL){ cta->set_L_if_correct_n(L,n); }
  }



  void low_end_multiply_and_compress(const MPS & other, RankCompressor &compressor, int sweep_start=0){
    if(a.size()<other.a.size()){
      std::cerr<<"MPS::multiply_front: a.size()<other.a.size()!"<<std::endl;
      exit(1);
    }
    for(size_t i=0; i<other.a.size(); i++){
std::cout<<"forward sweep "<<i<<"/"<<other.a.size()<<std::endl;

      a[i].multiply(other.a[i]);

      if(i>sweep_start){
        bool skip=false; 
        if(a[i-1].dim_d1==1 && a[i-1].dim_d2==1)skip=true;
        if(other.a[i-1].dim_d1==1 && other.a[i-1].dim_d2==1)skip=true;
        sweep_block_low_to_high(i-1, compressor);
      }
    }
    for(int i=other.a.size()-1; i>1; i--){
std::cout<<"backward sweep "<<i<<"/"<<other.a.size()<<std::endl;
      sweep_block_high_to_low(i, compressor);
    }
  }


  void high_end_multiply_and_compress(const MPS & other, RankCompressor &compressor, double keep_weight=0., Compress_Trafo_At *cta=NULL){
    if(a.size()<other.a.size()){
      std::cerr<<"MPS::high_end_multiply_and_compress: a.size()<other.a.size()!"<<std::endl;
      exit(1);
    }

    int start=((int)a.size())-((int)other.a.size());
    
    for(int i=((int)other.a.size())-1; i>=0; i--){
      a[i+start].multiply(other.a[i]);

      if(i<((int)other.a.size())-1){ 
        if(cta!=NULL && i+start+1==cta->n+1){
          sweep_block_high_to_low(i+start+1, compressor, keep_weight, &cta->L);
        }else{
          sweep_block_high_to_low(i+start+1, compressor, keep_weight);
        } 
      }
    }
    for(size_t i=start; i<a.size()-1; i++){
      if(cta!=NULL && i==cta->n){
        sweep_block_low_to_high(i, compressor, keep_weight, &cta->R);
      }else{
        sweep_block_low_to_high(i, compressor, keep_weight);
      }
    }
  }


  void multiply_and_compress(const MPS & other, RankCompressor &compressor){
    if(a.size()!=other.a.size()){
      std::cerr<<"MPS::multiply_front: a.size()<other.a.size()!"<<std::endl;
      exit(1);
    }
    low_end_multiply_and_compress(other, compressor);
  }

  void low_end_multiply(const MPS & other){
    int sz;
    if(a.size()<other.a.size()){
      sz=a.size();
    }else{
      sz=other.a.size();
    }
    for(size_t i=0; i<sz; i++){
      a[i].multiply(other.a[i]);
    }
  }
  virtual void sweep(RankCompressor &compressor, double keep_weight=0.){
     for(int n=0; n<(int)a.size()-1; n++){
       sweep_block_low_to_high(n, compressor, keep_weight);
     }
     for(int n=(int)a.size()-1; n>1; n--){
       sweep_block_high_to_low(n, compressor, keep_weight);
     }
  }

  void insert(int pos, const MPS_Matrix &b, int n_insert=1){
    int oldsize=a.size();
    std::vector<MPS_Matrix> mvec(oldsize+n_insert);
    for(size_t n=0; n<pos; n++){
      mvec[n].swap(a[n]);
    }
    for(size_t n=pos; n<pos+n_insert; n++){
      mvec[n].copy(b);
    }
    for(size_t n=pos; n<oldsize; n++){
      mvec[n_insert+n].swap(a[n]);
    }
    
    a.clear();
    a.resize(oldsize+n_insert);
    for(size_t n=0; n<oldsize+n_insert; n++){
      a[n].swap(mvec[n]);
    }
  }
  void write_binary(std::ostream &ofs)const{
    int sz=a.size();
std::cout<<"sz: "<<sz<<std::endl;
    ofs.write((char*)&sz, sizeof(int));
    for(int i=0; i<sz; i++){
      ofs.write((char*)&a[i].dim_i, sizeof(int));
      ofs.write((char*)&a[i].dim_d1, sizeof(int));
      ofs.write((char*)&a[i].dim_d2, sizeof(int));
    }
    for(int i=0; i<sz; i++){
      ofs.write((char*)a[i].mem, sizeof(std::complex<double>)*a[i].dim_i*a[i].dim_d1*a[i].dim_d2);
    }
  }
  void write_binary(const std::string &filename)const{
    std::ofstream ofs(filename.c_str(), std::ios::binary);
    write_binary(ofs);
  }
  virtual void read_binary(std::istream &ifs){
    int itmp;
    ifs.read((char*)&itmp, sizeof(int));
    std::cout<<"Reading size: "<<itmp<<std::endl;
    a.resize(itmp);
    for(int i=0; i<a.size(); i++){
      int dim_i, dim_d1, dim_d2;
      ifs.read((char*)&dim_i, sizeof(int));
      ifs.read((char*)&dim_d1, sizeof(int));
      ifs.read((char*)&dim_d2, sizeof(int));
      a[i].resize(dim_i, dim_d1, dim_d2);
    }
    for(int i=0; i<a.size(); i++){
      ifs.read((char*)a[i].mem, sizeof(std::complex<double>)*a[i].dim_i*a[i].dim_d1*a[i].dim_d2);
    }
  }
  virtual void read_binary(const std::string &filename){
    std::ifstream ifs(filename.c_str(), std::ios::binary);
    read_binary(ifs);
  }
  void copy(const MPS &other){
    a.clear();
    a.resize(other.a.size());
    for(size_t i=0; i<a.size(); i++){
      a[i].copy(other.a[i]);
    }
  }
  void swap(MPS &other){
    a.swap(other.a);
  }
  MPS &operator=(const MPS &other){
    copy(other);
    return *this;
  }
  MPS(const MPS &other){
    copy(other);
  }
  MPS(){}
};


#endif
