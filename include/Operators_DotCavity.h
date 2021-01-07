#ifndef OPERATORS_DOT_CAVITY_DEFINED_H
#define OPERATORS_DOT_CAVITY_DEFINED_H

#include "Operators.h"

class Operators_DotCavity{
public:
  int Nx_max;

  int get_dim()const{return 2*Nx_max+1;}
  //enumerate states: |G,0>, |X,0>, |G,1>, |X,1>, ...

  //with respect to dot state:
  Eigen::MatrixXcd ketbra(int i, int j){
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(get_dim(),get_dim());
    for(int nx=0; nx<Nx_max; nx++){ 
      mat(i+2*nx,j+2*nx)=1;
    }
    if(i<1 && j<1)mat(2*Nx_max+i,2*Nx_max+j)=1;
    return mat;
  }
  Eigen::MatrixXcd sigma_x(){
    return ketbra(0,1)+ketbra(1,0);
  }

  //with respect to dot state and cavity:
  Eigen::MatrixXcd ketbra(int i, int n, int j, int m){
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(get_dim(),get_dim());
    if(i+2*n>=get_dim() || j+2*m>=get_dim()){
      std::cerr<<"DotCavity::ketbra: i+2*n>=get_dim() || j+2*m>=get_dim()!"<<std::endl;
      exit(1);
    }
    mat(i+2*n,j+2*m)=1;
    return mat;
  }
  Eigen::MatrixXcd proj_n(int n){
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(get_dim(),get_dim());
    for(int i=0; i<2; i++){
      if(i+2*n<get_dim()){
        mat(i+2*n,i+2*n)=1;
      }
    }
    return mat;
  }
  Eigen::MatrixXcd n(){
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(get_dim(),get_dim());
    for(int n=0; n<Nx_max+1; n++){
      for(int i=0; i<2; i++){
        if(i+2*n<get_dim()){
          mat(i+2*n,i+2*n)=n;
        }
      }
    }
    return mat;
  }
  Eigen::MatrixXcd a(){
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(get_dim(),get_dim());
    for(int n=0; n<Nx_max+1; n++){
      for(int i=0; i<2; i++){
        if(i+2*n<get_dim() && i+2*(n-1)>=0){
          mat(i+2*(n-1),i+2*n)=sqrt((double)n);
        }
      }
    }
    return mat;
  }
  Eigen::MatrixXcd adagger(){
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(get_dim(),get_dim());
    for(int n=0; n<Nx_max+1; n++){
      for(int i=0; i<2; i++){
        if(i+2*(n+1)<get_dim()){
          mat(i+2*(n+1),i+2*n)=sqrt((double)n+1.);
        }
      }
    }
    return mat;
  }
  Eigen::MatrixXcd dot_cavity_coupling(){
    return ketbra(0,1)*adagger()+ketbra(1,0)*a();
  }

  Eigen::MatrixXcd zero(){
    return Eigen::MatrixXcd::Zero(get_dim(),get_dim());
  }

    

  Operators_DotCavity(int Nx_max_=1) : Nx_max(Nx_max_){
  }


};


#endif
