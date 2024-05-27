#ifndef SINGLE_BATH_MODE_DEFINED_H
#define SINGLE_BATH_MODE_DEFINED_H

#include <vector>
#include <Eigen/Dense>
#include "FreePropagator.hpp"

/**
Structure to store the interaction between a N-level system with a 
single bath mode

*/

//Note: indices: H( nu1*dim_xi+ xi1, nu2*dim_xi+xi2)

namespace ACE{

class SingleBathMode{
public:
  bool ready;

  int N; //dimension of N-level system
  int M; //dimension of bath mode

  Eigen::MatrixXcd H; //Hamiltonian
  //Time evolution in Liouville space: (vector indices: alpha, tilde{alpha}
  std::vector<std::vector<Eigen::MatrixXcd> > A;
 
  inline int get_N()const{return N;}
  inline int get_M()const{return M;}

  void check_dimensions();
  
  void calculateA_fullexp(double dt);

  void calculateA_factor1(double dt);
  
  inline void calculateA(double dt, int factor=0){
    if(factor==1){
      calculateA_factor1(dt);  
    }else{
      calculateA_fullexp(dt);
    }
  }

  inline SingleBathMode(int N_=0, int M_=0, const Eigen::MatrixXcd &H_=Eigen::MatrixXcd())
   : N(N_), M(M_), H(H_){
    check_dimensions();
    ready=false;
  }

};


}//namespace
#endif
