#ifndef ACE_REDUCED_LIOUVILLE_BASIS_BOSON_FB_DEFINED_H
#define ACE_REDUCED_LIOUVILLE_BASIS_BOSON_FB_DEFINED_H

#include "otimes.hpp"
#include "Operators_Boson.hpp"
#include "ReducedLiouvilleBasis.hpp"

/**  Reduce to the following basis:
     - start with initial state vector (initial mode density matrix)
     - construct new vectors by applying b^dagger and b at most 'nr_climb' times

HERE: Forward-backward version: 
     Add b/b^dagger once to forward and once backward direction but then also 
     include application of b/b^dagger at both, forward and backward direction.
     This matches more closely the effect of application of H_JC on the env.


*/
namespace ACE{

class ReducedLiouvilleBasis_Boson_FB: public ReducedLiouvilleBasis{
public:
 
  void setup_from_initial(const Eigen::MatrixXcd &initial, int nr_climb, double threshold=1e-12);

  inline ReducedLiouvilleBasis_Boson_FB(const Eigen::MatrixXcd &initial, int nr_climb, double threshold=1e-12){
    setup_from_initial(initial, nr_climb, threshold);
  }
  inline ReducedLiouvilleBasis_Boson_FB(const Eigen::MatrixXcd &U_)
    : ReducedLiouvilleBasis(U_){
  }
  inline ReducedLiouvilleBasis_Boson_FB(){
    use_reduce=false;
  }
};


}//namespace
#endif
