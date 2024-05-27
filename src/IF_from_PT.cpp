#include "IF_from_PT.hpp"
#include "ProcessTensorElement.hpp"
#include "ProcessTensor.hpp"
#include "InfluenceFunctional_OD.hpp"

namespace ACE{

InfluenceFunctional_OD IF_from_PT(const ProcessTensor &PT){
 
  if(PT.size()<1){
    std::cerr<<"IF_from_PT: PT.size()<1!"<<std::endl;
    exit(1);
  }

  InfluenceFunctional_OD IF;
  IF.dict=PT[0].accessor.dict;

  for(int n=1; n<PT.size(); n++){

    IF.dict.join(PT[n].accessor.dict);
  }

  IF.tgrid.set_default(PT.size());
  IF.a.resize(PT.size());
  IF.c.resize(PT.size());
  IF.env_ops.resize(PT.size());
  for(int n=0; n<PT.size(); n++){
    IF.a[n]=PT[n].M;
    if(IF.dict!=PT[n].accessor.dict){
      PT[n].accessor.dict.expand_MPS_Matrix(IF.a[n]);
      IF.dict.reduce_MPS_Matrix(IF.a[n]);
    }
    IF.c[n]=PT[n].closure;
    IF.env_ops[n]=PT[n].env_ops;
  }
   
  return IF;
}

}//namespace
