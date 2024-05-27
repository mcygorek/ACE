#ifndef ACE_EDM_ONE_BODY_OBSERVABLE_DEFINED_H
#define ACE_EDM_ONE_BODY_OBSERVABLE_DEFINED_H

#include "EDM_State.hpp"

namespace ACE {

class EDM_OneBodyObservable{
public:
  Eigen::VectorXcd o;
  
  std::complex<double> apply(Eigen::VectorXcd &rho_reduced) const;

  
 
};

}//namespace
#endif
