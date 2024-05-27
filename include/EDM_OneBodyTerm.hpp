#ifndef ACE_EDM_ONE_BODY_TERM_DEFINED_H
#define ACE_EDM_ONE_BODY_TERM_DEFINED_H

#include "EDM_State.hpp"

namespace ACE {

class EDM_OneBodyTerm{
public:
  virtual void update(double t, double dt) = 0;
  virtual EDM_State apply_filtered(const EDM_State &in, int r, const EDM_Filter &filter) = 0;
};

}
#endif
