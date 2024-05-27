#ifndef ACE_EFFECTIVE_PROPAGATOR_DEFINED_H
#define ACE_EFFECTIVE_PROPAGATOR_DEFINED_H

#include "Parameters.hpp"
#include "TimeGrid.hpp"
#include "FreePropagator.hpp"
#include "ProcessTensorForwardList.hpp"
#include "DummyException.hpp"

namespace ACE{
 
class EffectivePropagator: public Propagator{


static Eigen::MatrixXcd calculate(int at_n, 
           Propagator &prop, ProcessTensorForwardList &PT,
           const Eigen::MatrixXcd & initial_rho, const TimeGrid &tgrid);


};//class
}//namespace
#endif
