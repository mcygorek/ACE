#ifndef ACE_MPG_SELECTOR_DEFINED_H
#define ACE_MPG_SELECTOR_DEFINED_H

#include <vector>
#include <memory>
#include "ModePropagatorGenerator.hpp"

namespace ACE{
class Parameters;

class MPG_Selector{
public:
  std::vector<std::shared_ptr<ModePropagatorGenerator> > mpgs;

  inline operator std::vector<std::shared_ptr<ModePropagatorGenerator> >()const{
    return mpgs;
  }

  void setup(Parameters &param);

  inline MPG_Selector(Parameters &param){
    setup(param);
  }
};

}//namespace
#endif
