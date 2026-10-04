#ifndef POWERURV_DEFINED_H_
#define POWERURV_DEFINED_H_

#include "PCH.hpp"

namespace ACE{

void PowerURV_tall(Eigen::MatrixXcd &A, Eigen::MatrixXcd &L, Eigen::MatrixXcd &R, double threshold, int debug_mode=0);

void PowerURV_flat(Eigen::MatrixXcd &A, Eigen::MatrixXcd &L, Eigen::MatrixXcd &R, double threshold, int debug_mode=0);

void PowerURV(Eigen::MatrixXcd &A, Eigen::MatrixXcd &L, Eigen::MatrixXcd &R, double threshold, int debug_mode=0);

}//namespace
#endif
