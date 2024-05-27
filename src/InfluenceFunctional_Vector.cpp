#include "InfluenceFunctional_Vector.hpp"
#include <Eigen/Core>
#include "SpectralDensity.hpp"
#include "MPS.hpp"
#include "DiagBB.hpp"
#include <fstream>

namespace ACE{

  void InfluenceFunctional_Vector::calculate(int n_max){
    if(n_max<0){
      std::cerr<<"n_max must not be negative!"<<std::endl;
      exit(1);
    }
    
    b.resize(n_max+1);
    for(int n=0; n<n_max+1; n++){
//      b[n]=diagBB.calculate_expS(n_max-n, get_dt());
      b[n]=diagBB.calculate_expS(n, get_dt());
    }   
  }


  void InfluenceFunctional_Vector::print(const std::string &fname)const{
    std::cerr<<"InfluenceFunctional_Vector: print not implemented yet!"<<std::endl;
    exit(1);
  }

}//namespace
