#ifndef RRQR_DEFINED_H
#define RRQR_DEFINED_H


#include <iosfwd>
#include <vector>
#include "Eigen_fwd.hpp"
//#include <Eigen/Core>
//#include <Eigen/SVD>
//#include <iostream>
//#include <cstdlib>
//#include <vector>
//#include <fstream>

namespace ACE{

void HouseholderMultiply(const Eigen::VectorXcd &v,Eigen::MatrixXcd &A, int start_col=0);

std::vector<int> get_reverse_pivot(const std::vector<int> &pivot, bool print_pivots=false);

/** performs a rank-reducing factorization of a tall matrix A into L and R
    A is modified in the process */
void QRP(Eigen::MatrixXcd &A, Eigen::MatrixXcd &Q, Eigen::MatrixXcd &R, std::vector<int> &pivot, double threshold, int debug_mode=0);


void QLV(Eigen::MatrixXcd &A, Eigen::MatrixXcd &Q, Eigen::MatrixXcd &L, Eigen::MatrixXcd &V,double threshold, int debug_mode=0);


void QLV_with_debug(Eigen::MatrixXcd &A, Eigen::MatrixXcd &Q, Eigen::MatrixXcd &L, Eigen::MatrixXcd &V,double threshold, int debug_mode=0, int times=1);


void RRQR(Eigen::MatrixXcd &A, Eigen::MatrixXcd &L, Eigen::MatrixXcd &R, double threshold, int debug_mode=0, int times=1);

}//namespace
#endif
