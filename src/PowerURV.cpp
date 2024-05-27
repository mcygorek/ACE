#include "PowerURV.hpp"
#include <Eigen/Core>
#include <Eigen/SVD>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <vector>

namespace ACE{

void PowerURV_tall(Eigen::MatrixXcd &A, Eigen::MatrixXcd &L, Eigen::MatrixXcd &R, double threshold, int debug_mode){

  if(A.rows()<A.cols()){
    std::cerr<<"PowerURV: Matrix A not tall!"<<std::endl;
    exit(1);
  }
  if(A.cols()<1){
    L.resize(A.rows(), A.cols());
    R.resize(A.cols(), A.cols());
    return;
  }
  if(A.cols()==1){
//    double n=A.col(0).norm();
//    if(n<1e-16)n=1.;
    double n=1.;
    L=A/n;
    R.resize(1,1); R(0,0)=n;
    return; 
  }

//  bool silent=!(debug_mode&1);
  bool test_result=(debug_mode&2);
//  bool print_pivots=(debug_mode&4);
  bool print_matrices=(debug_mode&8);
//  bool nopivot=(debug_mode &16);
//  double compare_zero=1e-16;

//  int threshold_reached=-1;
//  double threshold_remainder=0.;
  
  if(print_matrices){
    std::cout<<"Original matrix:"<<std::endl;
    std::cout<<A<<std::endl<<std::endl;
  }
  Eigen::MatrixXcd A_backup;
  if(test_result){
    A_backup=A;
  } 
 
  //----------------------


  Eigen::MatrixXcd G=Eigen::MatrixXcd::Random(A.cols(),A.cols()); 
  Eigen::MatrixXcd AA=A.adjoint()*A;
//  AA*=AA;

  Eigen::HouseholderQR<Eigen::MatrixXcd> qr_start(AA*G); 
  Eigen::MatrixXcd V=qr_start.householderQ();

  Eigen::HouseholderQR<Eigen::MatrixXcd> qr(A*V); 
  Eigen::MatrixXcd Q=qr.householderQ();
  Eigen::MatrixXcd Rtmp=Q.adjoint()*A*V;
 

  int k=Rtmp.cols();
  for(int i=0; i<Rtmp.cols(); i++){
    if(std::abs(Rtmp(i,i))<threshold){
      k=i;
      break;
    }
  }
  if(k<1){
    std::cerr<<"EigenQR: rank < 1"<<std::endl;
    exit(1);
  }
/*
  std::cout<<"rank: "<<k<<std::endl;

  std::cout<<"Rtmp diags: ";
  for(int i=0; i<Rtmp.cols(); i++){
    std::cout<<Rtmp(i,i)<<" ";
  }std::cout<<std::endl;
*/

/*   
  L=Q;
  R=Rtmp*V.adjoint();
*/
  L=Q.block(0,0,A.rows(),k);
  R=Rtmp.block(0,0,k,A.cols())*V.adjoint();
  

  

  //----------------------

  if(print_matrices){
    std::cout<<"L^dagger*L:"<<std::endl<<L.adjoint()*L<<std::endl<<std::endl;
    std::cout<<"L*R:"<<std::endl<<L*R<<std::endl<<std::endl;
  }

  if(test_result){
    Eigen::MatrixXcd B=L*R;
    double maxdiff=0.;
    for(int i=0; i<A.rows(); i++){
      for(int j=0; j<A.cols(); j++){
        double a=abs(A_backup(i,j)-B(i,j));
        if(a>maxdiff)maxdiff=a;
      }
    }
    std::cout<<"Maximal difference between A and QR: "<<maxdiff<<std::endl<<std::endl;
    if(maxdiff>threshold*100){
      std::ofstream ofs("DUMP_A.dat");
      for(int i=0; i<A.rows(); i++){
        for(int j=0; j<A.cols(); j++){
          ofs<<A(i,j)<<" ";
        }
        ofs<<std::endl;
      }
    }
  }
}

void PowerURV_flat(Eigen::MatrixXcd &A, Eigen::MatrixXcd &L, Eigen::MatrixXcd &R, double threshold, int debug_mode){
  Eigen::MatrixXcd Ltmp, Rtmp;
  Eigen::MatrixXcd Atmp=A.adjoint();

  PowerURV_tall(Atmp, Ltmp, Rtmp, threshold, debug_mode);
  R=Ltmp.adjoint();
  L=Rtmp.adjoint();
}

void PowerURV(Eigen::MatrixXcd &A, Eigen::MatrixXcd &L, Eigen::MatrixXcd &R, double threshold, int debug_mode){
  if(A.rows()>=A.cols()){
    PowerURV_tall(A, L, R, threshold, debug_mode);
  }else{
    PowerURV_flat(A, L, R, threshold, debug_mode);
  }
}

}//namespace
