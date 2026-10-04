#include "CombElement.hpp"
#include "DummyException.hpp"
#include <sstream>
#include <iostream>

namespace ACE{

void CombElement::sweep_d2_d4(const TruncatedSVD &trunc, Eigen::MatrixXcd &P, bool is_last, int verbosity){

    Eigen::MatrixXcd A;
    if(P.rows()!=0){
      if(P.rows()!=f.dim_d2){
        std::cerr<<"CombElement::sweep_d2_d4: P.rows()!=f.dim_d2!"<<std::endl;
        throw DummyException();
      }
      f.multiply_d2(P);
    }
    if(is_last)return;

    A = f.get_Matrix_d1d2d3_d4();
    TruncatedSVD_Return ret = trunc.compress(A);
//    std::cout<<"["<<l<<"]: singular values: "<<ret.sigma.transpose()<<std::endl;
    double keep=ret.sigma(0);
    if(trunc.keep!=0.)keep=trunc.keep;
#ifdef PRINT_KEEP
std::cout<<"keep: "<<keep<<std::endl;
#endif
    f.set_from_Matrix_d1d2d3_d4( ret.U*keep, f.dim_d2, f.dim_d3);
    P=((ret.sigma/keep).asDiagonal()*ret.Vdagger).transpose();
    if(verbosity>0){std::cout<<"forward sweep: "<<A.rows()<<","<<A.cols()<<" -> "<<P.cols()<<std::endl;}
}

void CombElement::sweep_d4_d2(const TruncatedSVD &trunc, Eigen::MatrixXcd &P, bool is_last, int verbosity){

    Eigen::MatrixXcd A;
    if(P.rows()!=0){
      if(P.rows()!=f.dim_d4){
        std::cerr<<"CombElement::sweep_d4_d2: P.rows()!=f.dim_d4!"<<std::endl;
        throw DummyException();
      }
      f.multiply_d4(P);
    }
    if(is_last)return;

    A = f.get_Matrix_d2_d1d3d4();
    TruncatedSVD_Return ret = trunc.compress(A);
//    std::cout<<"["<<l<<"]: singular values: "<<ret.sigma.transpose()<<std::endl;
    double keep=ret.sigma(0);
    if(trunc.keep!=0.)keep=trunc.keep;
#ifdef PRINT_KEEP
std::cout<<"keep: "<<keep<<std::endl;
#endif
    f.set_from_Matrix_d2_d1d3d4( ret.Vdagger*keep, f.dim_d3, f.dim_d4);
    P=(ret.U*(ret.sigma/keep).asDiagonal());
    if(verbosity>0){std::cout<<"backward sweep: "<<A.rows()<<","<<A.cols()<<" -> "<<P.cols()<<std::endl;}
}


void CombElement::read_binary(std::istream &is){
  f.read_binary(is);
  closure_d1=binary_read_EigenMatrixXcd(is,"closure_d1");

  { int sz=binary_read_int(is, "tooth exists");
    if(sz<0){ tooth.reset();} 
    else{
      tooth=std::make_unique<ProcessTensorBuffer>();
      tooth->resize(sz);
      for(int i=0; i<tooth->n_tot; i++){  
        tooth->get(i).read_binary(is);
      }
    }
  }

//  { int sz=binary_read_int(is, "b.size()");
//    b.resize(sz); }
//  for(size_t i=0; i<b.size(); i++){
//    b[i]=binary_read_EigenMatrixXcd(is,"b elements");
//  }

}
void CombElement::write_binary(std::ostream &os)const{
  f.write_binary(os);
  binary_write_EigenMatrixXcd(os, closure_d1);

  if(!tooth){ binary_write_int(os, -1); }
  else{ 
    binary_write_int(os, tooth->n_tot); 
    for(int i=0; i<tooth->n_tot; i++){  
      tooth->get(i).write_binary(os);
    }
  }

//  binary_write_int(os, b.size());
//  for(size_t i=0; i<b.size(); i++){
//    binary_write_EigenMatrixXcd(os,b[i]);
//  }

} 
CombElement::~CombElement(){};
}//namespace
