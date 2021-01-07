#ifndef FULL_INFLUENCE_FUNCTIONAL_DENSE_DEFINED_H
#define FULL_INFLUENCE_FUNCTIONAL_DENSE_DEFINED_H

#include "Tensor.h"
#include "InfluenceFunctional.h"

class FullInfluenceFunctional_Dense: public Tensor_Dense{
public:


  FullInfluenceFunctional_Dense(const InfluenceFunctional &IF){
 
    int NL=IF.get_dim()*IF.get_dim();
    int n_max=IF.get_n_max();
    resize(n_max+1,NL);

    int bs=1; for(int i=0; i<n_max+1; i++)bs*=NL;
    for(int b=0; b<bs; b++){
      mem[b]=IF.ten[n_max][b];
    }
    for(int n=0; n<n_max;n++){
      int bs1=1; for(int i=0; i<n+1; i++)bs1*=NL;
      int bs2=bs/bs1;
      for(int b=0; b<bs1; b++){
        for(int b2=0; b2<bs2; b2++){
          mem[b2*bs1+b]*=IF.ten[n][b];
        }
      }
    }
  }
  FullInfluenceFunctional_Dense(){
  }
};


#endif
