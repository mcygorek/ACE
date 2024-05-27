#ifndef ACE_EDM_ONE_BODY_TERM_FPROP_DEFINED_H
#define ACE_EDM_ONE_BODY_TERM_FPROP_DEFINED_H

#include "EDM_OneBodyTerm.hpp"
#include "FreePropagator.hpp"
#include "Parameters.hpp"
#include <memory>

namespace ACE {

class EDM_OneBodyTerm_FreePropagator: public EDM_OneBodyTerm{
public:
  std::shared_ptr<FreePropagator> fprop;

  inline operator bool()const{
    return (bool)fprop;
  }

  virtual void update(double t, double dt);
  virtual EDM_State apply_filtered(const EDM_State &in, int r, const EDM_Filter &filter);
 
  void setup(Parameters &param, int site);
 
  EDM_OneBodyTerm_FreePropagator(){}
  EDM_OneBodyTerm_FreePropagator(Parameters &param, int site){
    setup(param, site);
  }
};

}//namespace
#endif
