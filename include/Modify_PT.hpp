#ifndef ACE_MODIFY_PT_DEFINED_H
#define ACE_MODIFY_PT_DEFINED_H

#include "InfluenceFunctional_OD.hpp"

/* Contains routines to modify existing process tensors, e.g.,
   Coarse graining or multiplication with a system propagator.
 
   Operations are usually in-place. 
*/

namespace ACE{
namespace Modify_PT{

// PT with time step dt -> PT with time step (n_coarse*dt)
void coarse_grain(InfluenceFunctional_OD &IF, int n_coarse, double dict_zero=1e-20);

void apply_system_propagator(InfluenceFunctional_OD &IF, Propagator &prop, double ta, double dt, double dict_zero=1e-20);
   
}//namespace
}//namespace


#endif
