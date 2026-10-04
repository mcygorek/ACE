#ifndef ACE_FENNA_MATTHEWS_OLSON_DEFINED_H_
#define ACE_FENNA_MATTHEWS_OLSON_DEFINED_H_

#include <Eigen/Dense>
#include "Constants.hpp"


namespace ACE{
//see: Adolphs & Renger, Biophysical Journal 91, 2778 (2006)

inline Eigen::MatrixXcd H_FMO_monomer(){
  Eigen::MatrixXcd H(7,7);
  H <<   240., -87.7, 5.5, -5.9, 6.7, -13.7, -9.9,
         -87.7, 315., 30.8, 8.2, 0.7, 11.8, 4.3,
         5.5, 30.8, 0., -53.5, -2.2, -9.6, 6.0,
         -5.9, 8.2, -53.5, 130, -70.7, -17.0, -63.3,
         6.7, 0.7, -2.2, -70.7, 285., 81.1, -1.3,
         -13.7, 11.8, -9.6, -17.0, 81.1, 435., 39.7,
         -9.9, 4.3, 6.0, -63.3, -1.3, 39.7, 245;
  return H*inv_cm_in_meV;
}
inline double H_FMO_monomer_shift(){ return 12205*inv_cm_in_meV; }

inline Eigen::MatrixXcd H_FMO_trimer(){
  Eigen::MatrixXcd H(7,7);
  H <<   200., -87.7, 5.5, -5.9, 6.7, -13.7, -9.9,
         -87.7, 320., 30.8, 8.2, 0.7, 11.8, 4.3,
         5.5, 30.8, 0., -53.5, -2.2, -9.6, 6.0,
         -5.9, 8.2, -53.5, 110., -70.7, -17.0, -63.3,
         6.7, 0.7, -2.2, -70.7, 270., 81.1, -1.3,
         -13.7, 11.8, -9.6, -17.0, 81.1, 420., 39.7,
         -9.9, 4.3, 6.0, -63.3, -1.3, 39.7, 230;
  return H*inv_cm_in_meV;
}
inline double H_FMO_trimer_shift(){ return 12210*inv_cm_in_meV; }

}//namespace
#endif
