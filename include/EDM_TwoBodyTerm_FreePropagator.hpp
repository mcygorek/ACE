#ifndef ACE_EDM_TWO_BODY_TERM_FPROP_DEFINED_H
#define ACE_EDM_TWO_BODY_TERM_FPROP_DEFINED_H

#include "EDM_TwoBodyTerm.hpp"
#include "FreePropagator.hpp"
#include "Parameters.hpp"

namespace ACE {
//          ------
//alpha01---|    |---alpha00
//          |    |
//alpha11---|    |---alpha10
//          ------

class EDM_TwoBodyTerm_FreePropagator: public EDM_TwoBodyTerm{
public:
  std::shared_ptr<FreePropagator> fprop;

  inline operator bool()const{
    return (bool)fprop;
  }

  virtual void update(int n, const TimeGrid &tgrid);
  virtual EDM_State apply_filtered(const EDM_State &in, const std::pair<int,int> &r, const EDM_Filter &filter);

  void setup(Parameters &param, const std::pair<int,int> &site);

  EDM_TwoBodyTerm_FreePropagator(){}
  EDM_TwoBodyTerm_FreePropagator(Parameters &param, const std::pair<int,int> &site){
    setup(param, site);
  }

};
}
#endif
