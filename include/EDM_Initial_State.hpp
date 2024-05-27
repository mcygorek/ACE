#ifndef ACE_EDM_INITIAL_STATE_DEFINED_H
#define ACE_EDM_INITIAL_STATE_DEFINED_H

#include "EDM_State.hpp"
#include "Parameters.hpp"
#include <Eigen/Core>

namespace ACE {

class EDM_Initial_State{
public:
  EDM_State rhos;

  inline operator EDM_State&(){
    return rhos;
  }

  void set_default_closure(int r);

  void setup(Parameters &param);

  EDM_Initial_State(){}
  EDM_Initial_State(Parameters &param){
    setup(param);
  }

};
}
#endif
