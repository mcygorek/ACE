#ifndef HERMITIAN_LIOUVILLE_BASIS_DEFINED_H
#define HERMITIAN_LIOUVILLE_BASIS_DEFINED_H

#include <Eigen/Dense>

/* Liouville space basis, in the follwing order:

- Trace,
- diagonals (ortho. wrt. trace; omit ground state)
- hermitian combination of off-diagonals
- anti-hermitian combination of off-diagonals

*/

namespace ACE{

class HermitianLiouvilleBasis{
private:
  int N;
public:
  inline int get_N()const{return N;}
  inline void set_N(int N_){N=N_;}
  
  Eigen::VectorXcd operator() (int i) const;
  
  Eigen::MatrixXcd get_Matrix()const;

  inline HermitianLiouvilleBasis(int N_):N(N_){}
  inline virtual ~HermitianLiouvilleBasis(){}
};

}//namespace
#endif
