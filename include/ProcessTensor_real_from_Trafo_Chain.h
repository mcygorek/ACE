#ifndef ACE_PTREAL_FROM_TRAFO_CHAIN_DEFINED_H
#define ACE_PTREAL_FROM_TRAFO_CHAIN_DEFINED_H

#include "ProcessTensor_real.hpp"

namespace ACE{

class ProcessTensor_real_from_Trafo_Chain{
public:

  ProcessTensor_real PT;

  ProcessTensor_real get_PT(int nmax){
    if(PT.a.size()<3){
      std::cerr<<"ProcessTensor_real_from_Trafo_Chain: PT.a.size()<3!"<<std::endl;
      exit(1);
    }
    ProcessTensor_real PT2;
    PT2.dict=PT.dict;
    PT2.HLU=PT.HLU;
    PT2.c=std::vector<Eigen::VectorXcd>(n_max, PT.c[0]);
    PT2.env_ops=std::vector<std::vector<Eigen::VectorXcd> >(n_max, PT.env_ops[0]);
    PT2.a=std::vector<MPS_Matrix_real>(n_max, PT.a[1]);
    PT2.a[0]=PT.a[0];
    PT2.a.back()=PT.a[2];

    return PT2;
  }

  void calculate(ModePropagatorGenerator &mpg, Trafo_Chain &trafo_chain){
  }

  ProcessTensor_real_from_Trafo_Chain(ModePropagatorGenerator &mpg, Trafo_Chain &trafo_chain){
    calculate(mpg, trafo_chain);
  }
};

}//namespace
#endif
