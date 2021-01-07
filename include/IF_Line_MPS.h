#ifndef IF_LINE_MPS_DEFINED_H
#define IF_LINE_MPS_DEFINED_H

#include <Eigen/Core>
#include "DiagBB.h"
#include "MPS.h"

class IF_Line_MPS: public MPS{
public:
  
  void calculate(int n, double dt, DiagBB & diagBB, int n_tot=0){
    if(n_tot==0)n_tot=n;
    int NL=diagBB.get_dim()*diagBB.get_dim();
   
    MPS::resize(n_tot+1, NL);
    if(n<0){
      std::cerr<<"IF_Line_MPS::calculate: n<0!"<<std::endl;
      exit(1);
    }

    int n_one=n_tot-n;
//    std::cout<<"n_one: "<<n_one<<std::endl;
    for(int i=0; i<n_one; i++){
      a[i].resize(NL,1,1);
      a[i].fill(1.);
    }
    if(n==0){  //special case: matrix dimensions NL,1,1
      Eigen::MatrixXcd expS=diagBB.calculate_expS(0,dt);
      for(int i=0; i<NL; i++){
        a[n_one](i,0,0)=expS(i,i);
      }
    }else{
      a[n_one].resize(NL,1,NL);
      a[n_one].set_zero();
      Eigen::MatrixXcd expS=diagBB.calculate_expS(n,dt);
      for(int i=0; i<NL; i++){
        for(int j=0; j<NL; j++){
          a[n_one](i, 0, j)=expS(i,j);
        }
      }

      for(int k=1; k<n; k++){   
        a[n_one+k].resize(NL,NL,NL);
        a[n_one+k].set_zero();
        expS=diagBB.calculate_expS(n-k,dt);
        for(int i=0; i<NL; i++){
          for(int j=0; j<NL; j++){
            a[n_one+k](i, j, j)=expS(i,j);
          }
        } 
      }


      a[n_tot].resize(NL,NL,1);
      a[n_tot].set_zero();
      expS=diagBB.calculate_expS(0,dt);
      for(int i=0; i<NL; i++){
          a[n_tot](i, i, 0)=expS(i,i);
      }
    }
/*
    std::cout<<"_"<<std::endl;
    for(int i=0; i<n_tot+1;i++){
      std::cout<<a[i].dim_d1<<" ";
    }
*/
    std::cout<<std::endl;
  }

  IF_Line_MPS(int n, double dt, DiagBB & diagBB, int n_tot=0){
    calculate(n, dt, diagBB, n_tot);
  }

  IF_Line_MPS(){
  }
  
};


#endif
