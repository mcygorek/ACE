#ifndef ACE_DIMENSION_EXTENDER_DEFINED_H
#define ACE_DIMENSION_EXTENDER_DEFINED_H

#include <iostream>
#include "otimes.hpp"

/**
  Extend matrix acting on a lower subspace to a higher subspace. 
  E.g., H_S -> (H_S otimes Id_2) 

  Dimensions defined in Hilbert space; extension of Hilbert or Liouville space.
*/

namespace ACE{

class DimensionExtender{
public:
  //conceptualize factorization of input matrix
  //e.g.: H_in = H0 x H1 x H2, where dim(H0)=3, dim(H1)=2, dim(H2)=4 -> {3,2,4}
  std::vector<int> dims_in; 

  //dimensions to extend: 
  // H_in = H0 x H1 x H2 -> E0 x H0 x E1 x H1 x E2 x H2 x E3
  // where Ei=Id_{dim_ext[i]}   => dims_ext.size()=dims_in.size()+1;
  std::vector<int> dims_ext; 

  
  int dim_in()const{
    int indim=1;
    for(size_t i=0; i<dims_in.size(); i++)indim*=dims_in[i];
    return indim;
  }
  int dim_ext()const{
    int extdim=1;
    for(size_t i=0; i<dims_ext.size(); i++)extdim*=dims_ext[i];
    return extdim;
  }
  int dim_out()const{
    return dim_in()*dim_ext();
  }
  void check_consistency(){
    if(dims_in.size()<1){
      std::cerr<<"Error: DimensionExtender: dims_in.size()<1!"<<std::endl;
      exit(1);
    }
    if(dims_ext.size()!=dims_in.size()+1){
      std::cerr<<"Error: DimensionExtender: dims_ext.size()!=dims_in.size()+1!"<<std::endl;
      exit(1);
    }
    for(size_t i=0; i<dims_in.size(); i++){
      if(dims_in[i]<1){
        std::cerr<<"Error: DimensionExtender: dims_in["<<i<<"]<1!"<<std::endl;
        exit(1);
      }
    }
    for(size_t i=0; i<dims_ext.size(); i++){
      if(dims_ext[i]<1){
        std::cerr<<"Error: DimensionExtender: dims_ext["<<i<<"]<1!"<<std::endl;
        exit(1);
      }
    }
  }
  void check_consistency_Hilbert(const Eigen::MatrixXcd &M){
    check_consistency();
    if(M.rows()!=M.cols()){
      std::cerr<<"Error: DimensionExtender: M.row()!=M.col()!"<<std::endl;
      exit(1);
    }

    int indim=dim_in();
    if(indim!=M.rows()){
      std::cerr<<"Error: DimensionExtender: indim!=M.row() ("<<indim<<" vs. "<<M.rows()<<"!"<<std::endl;
      exit(1);
    }
  }
  void check_consistency_Liouville(const Eigen::MatrixXcd &M){
    check_consistency();
    if(M.rows()!=M.cols()){
      std::cerr<<"Error: DimensionExtender: M.row()!=M.col()!"<<std::endl;
      exit(1);
    }

    int indim=dim_in();
    int NH=sqrt(M.rows());
    if(NH*NH!=M.rows()){
      std::cerr<<"Error: DimensionExtender: NH*NH!=M.rows() ("<<NH*NH<<" vs. "<<M.rows()<<")!"<<std::endl;
      exit(1);
    }
    if(indim!=NH){
      std::cerr<<"Error: DimensionExtender: dim_in() != NH ("<<indim<<" vs. "<<M.rows()<<")!"<<std::endl;
      exit(1);
    }
  }

 
  //get big inner index and ext index from output index (last arg. optional to speed up)
  std::pair<int, int> in_ext_index(int iout, int dout=0){
    int iin=0; 
    int iext=0;
    if(dout<1)dout=dim_out();
    for(size_t i=0; i<dims_in.size(); i++){
      iin*=dims_in[i]; 
      iext*=dims_ext[i];

      dout=dout/dims_ext[i];
      iext+=iout/dout;
      iout=iout%dout;
  
      dout=dout/dims_in[i];
      iin+=iout/dout;
      iout=iout%dout;
    } 
    iext*=dims_ext.back();
    iext+=iout;
    return std::make_pair(iin, iext);
  }
  int in_index(int iout, int dout=0){
    return in_ext_index(iout, dout).first;
  }

  // gets inner indices of two outer indices; returns 'false' if extendend parts don't match
  bool in_index_delta(int iout1, int iout2, int &iin1, int &iin2, int dout=0){
    iin1=0; iin2=0;
    if(dout<1)dout=dim_out();
    std::pair<int, int> in_ext1=in_ext_index(iout1, dout);
    std::pair<int, int> in_ext2=in_ext_index(iout2, dout);
    iin1=in_ext1.first;
    iin2=in_ext2.first;
    return in_ext1.second==in_ext2.second;
  }
/*
  bool in_index_delta(int iout1, int iout2, int &iin1, int &iin2, int dout=0){
    iin1=0; iin2=0;
    if(dout<1)dout=dim_out();

    if(iout1%dims_ext.back() != iout2%dims_ext.back())return false;
    for(size_t i=0; i<dims_in.size(); i++){
      if(i>0)iin1*=dims_in[i-1]; 
      if(i>0)iin2*=dims_in[i-1]; 
      dout=dout/dims_ext[i];
      if(iout1/dout!=iout2/dout)return false;

      iout1=iout1%dout;
      iout2=iout2%dout;
      dout=dout/dims_in[i];
      iin1+=iout1/dout;
      iin2+=iout2/dout;
    }
    return true;
  }
*/
 
  void initialize(){
    dims_in.resize(1,1);
    dims_ext.resize(2,1);
  }
  void initialize(int Nfront, int N1, int Nback){
    dims_in.resize(1, N1);
    dims_ext.resize(2); 
    dims_ext[0]=Nfront; 
    dims_ext[1]=Nback; 
  }
  void initialize(int Nfront, int N1, int Nmid, int N2, int Nback){
    dims_in.resize(2); dims_in[0]=N1; dims_in[1]=N2;
    dims_ext.resize(3); 
    dims_ext[0]=Nfront; 
    dims_ext[1]=Nmid; 
    dims_ext[2]=Nback; 
  }

  Eigen::MatrixXcd get_Hilbert(const Eigen::MatrixXcd &M){
    check_consistency_Hilbert(M);
    int dout=dim_out();

    Eigen::MatrixXcd A=Eigen::MatrixXcd::Zero(dout, dout);

    //strategy: loop through output dimensions; get corresponding input index
    int rin, cin;
    for(int r=0; r<dout; r++){
      for(int c=0; c<dout; c++){
        if(in_index_delta(r, c, rin, cin, dout)){
          A(r,c)=M(rin,cin);    
        }
      }
    }
    return A;
  }
  Eigen::MatrixXcd get_Liouville(const Eigen::MatrixXcd &M){
    check_consistency_Liouville(M);
    int din=dim_in();
    int dout=dim_out();

    Eigen::MatrixXcd A=Eigen::MatrixXcd::Zero(dout*dout, dout*dout);

    //strategy: loop through output dimensions; get corresponding input index

    int r1in, c1in;
    int r2in, c2in;
    for(int r1=0; r1<dout; r1++){
      for(int c1=0; c1<dout; c1++){
        if(in_index_delta(r1, c1, r1in, c1in, dout)){
          for(int r2=0; r2<dout; r2++){
            for(int c2=0; c2<dout; c2++){
              if(in_index_delta(r2, c2, r2in, c2in, dout)){
                A(r1*dout+r2,c1*dout+c2) = M(r1in*din+r2in,c1in*din+c2in);    
              }
            }
          }
        }
      }
    }

    return A;
  }
  DimensionExtender(){
    initialize();
  }
  DimensionExtender(int Nfront, int N1, int Nback){
    initialize(Nfront, N1, Nback);
  }
  DimensionExtender(int Nfront, int N1, int Nmid, int N2, int Nback){
    initialize(Nfront, N1, Nmid, N2, Nback);
  }
};

}//namespace
#endif
