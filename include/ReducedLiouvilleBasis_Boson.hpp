#ifndef ACE_REDUCED_LIOUVILLE_BASIS_BOSON_DEFINED_H
#define ACE_REDUCED_LIOUVILLE_BASIS_BOSON_DEFINED_H

#include "ReducedLiouvilleBasis.hpp"

/**  Reduce to the following basis:
     - start with initial state vector (initial mode density matrix)
     - construct new vectors by applying b^dagger and b at most 'nr_climb' times
*/
namespace ACE{

class ReducedLiouvilleBasis_Boson: public ReducedLiouvilleBasis{
public:
 
  void setup_from_initial(const Eigen::MatrixXcd &initial, int nr_climb, double threshold=1e-12);

  inline ReducedLiouvilleBasis_Boson(const Eigen::MatrixXcd &initial, int nr_climb, double threshold=1e-12){
    setup_from_initial(initial, nr_climb, threshold);
  }
  inline ReducedLiouvilleBasis_Boson(const Eigen::MatrixXcd &U_)
    : ReducedLiouvilleBasis(U_){
  }
  inline ReducedLiouvilleBasis_Boson(){
    use_reduce=false;
  }
};


}//namespace
#endif
