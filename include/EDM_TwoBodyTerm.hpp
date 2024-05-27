#ifndef ACE_EDM_TWO_BODY_TERM_DEFINED_H
#define ACE_EDM_TWO_BODY_TERM_DEFINED_H

#include "EDM_State.hpp"
#include "TimeGrid.hpp"

namespace ACE {
//          ------
//alpha01---|    |---alpha00
//          |    |
//alpha11---|    |---alpha10
//          ------

class EDM_TwoBodyTerm{
public:
  virtual void update(int n, const TimeGrid &tgrid) = 0;

  virtual EDM_State apply_filtered(const EDM_State &in, const std::pair<int,int> &r, const EDM_Filter &filter) = 0;
};
}
#endif
