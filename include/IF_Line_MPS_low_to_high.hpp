#ifndef IF_LINE_MPS_LOW_TO_HIGH_DEFINED_H
#define IF_LINE_MPS_LOW_TO_HIGH_DEFINED_H

#include "DiagBB.hpp"
#include "MPS.hpp"

namespace ACE{

class IF_Line_MPS_low_to_high: public MPS{
public:
  
  void calculate(int n, double dt, DiagBB & diagBB, int n_tot=0);

  IF_Line_MPS_low_to_high(int n, double dt, DiagBB & diagBB, int n_tot=0){
    calculate(n, dt, diagBB, n_tot);
  }

  IF_Line_MPS_low_to_high(){
  }
  
};

}//namespace
#endif
