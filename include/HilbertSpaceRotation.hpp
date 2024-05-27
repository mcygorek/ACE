#pragma once
#ifndef ACE_HILBERT_SPACE_ROTATION_DEFINED_H
#define ACE_HILBERT_SPACE_ROTATION_DEFINED_H

#include <Eigen/Dense>
/*
   Purpose: Store and apply rotations of the Hilbert space to propagators,
   initial values, and output operators
*/

namespace ACE{
class HilbertSpaceRotation{
public:

  bool use;  //do nothing, if false
  Eigen::MatrixXcd U; //Columns: Basis vectors

  inline bool used()const{return use;}

  Eigen::MatrixXcd apply(const Eigen::MatrixXcd & A)const;

  Eigen::MatrixXcd apply_Liouville(const Eigen::MatrixXcd & L)const;
  
  void setup_by_diagonalizing(const Eigen::MatrixXcd & A);
  
  inline HilbertSpaceRotation(){
    use=false;
  }

};

}//namespace
#endif
