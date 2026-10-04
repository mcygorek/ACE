#pragma once
#ifndef MPO_APPLY_OPERATOR_DEFINED_H
#define MPO_APPLY_OPERATOR_DEFINED_H

#include "MPO.hpp"
#include "Comb.hpp"
#include "Parameters.hpp"

namespace ACE{


struct MPO_ApplyOperators{
  std::vector<std::pair< int, MPO > > list;
  
  void apply(int step, Comb & comb, const TruncatedSVD &trunc, int verbosity=1);

  void setup(Parameters & param);

  MPO_ApplyOperators(Parameters & param){setup(param);}
  MPO_ApplyOperators(){}
  ~MPO_ApplyOperators(){}
};
}//namespace
#endif
