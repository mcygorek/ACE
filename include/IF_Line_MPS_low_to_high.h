#ifndef IF_LINE_MPS_LOW_TO_HIGH_DEFINED_H
#define IF_LINE_MPS_LOW_TO_HIGH_DEFINED_H

#include <Eigen/Core>
#include "DiagBB.h"
#include "MPS.h"

class IF_Line_MPS_low_to_high: public MPS{
public:
  
  void calculate(int n, double dt, DiagBB & diagBB, int n_tot=0){
    if(n_tot<n)n_tot=n;
    int NL=diagBB.get_dim()*diagBB.get_dim();
   
    MPS::resize(n_tot, NL);
    if(n<=0){
      std::cerr<<"IF_Line_MPS::calculate: n<=0!"<<std::endl;
      exit(1);
    }

    if(n_tot==1){  //special case: matrix dimensions NL,1,1
      a[0].resize(NL,1,1);
      Eigen::MatrixXcd expS=diagBB.calculate_expS(0,dt);
      for(int i=0; i<NL; i++){
        a[0](i,0,0)=expS(i,i);
      }
    }else{
      a[0].resize(NL,1,NL);
      a[0].set_zero();
      Eigen::MatrixXcd expS=diagBB.calculate_expS(0,dt);
      for(int i=0; i<NL; i++){
          a[0](i, 0, i)=expS(i,i);
      }

      for(int k=1; k<n-1; k++){   
        a[k].resize(NL,NL,NL);
        a[k].set_zero();
        expS=diagBB.calculate_expS(k,dt);
        for(int i=0; i<NL; i++){
          for(int j=0; j<NL; j++){
            a[k](i, j, j)=expS(i,j);
          }
        } 
      }


      a[n-1].resize(NL,NL,1);
      a[n-1].set_zero();
      expS=diagBB.calculate_expS(n-1,dt);
      for(int i=0; i<NL; i++){
        for(int j=0; j<NL; j++){
          a[n-1](i, j, 0)=expS(i,j);
        }
      }
    }

    for(int i=n; i<n_tot; i++){
      a[i].resize(NL,1,1);
      a[i].fill(1.);
    }
    std::cout<<std::endl;
  }

  IF_Line_MPS_low_to_high(int n, double dt, DiagBB & diagBB, int n_tot=0){
    calculate(n, dt, diagBB, n_tot);
  }

  IF_Line_MPS_low_to_high(){
  }
  
};


#endif
