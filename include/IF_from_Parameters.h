#ifndef IF_FROM_PARAMETERS_DEFINED_H
#define IF_FROM_PARAMETERS_DEFINED_H

#include "IF_OD_Plan.h"

std::vector<Smart_Ptr<InfluenceFunctional_OD> > IF_from_Parameters(Parameters &param, int check_sys_dim=0){
  IF_OD_Plan plan(param, check_sys_dim);
  return plan.execute();
}

#endif
