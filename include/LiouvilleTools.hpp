#pragma once
#ifndef ACE_LIOUVILLE_TOOLS_DEFINED_H
#define ACE_LIOUVILLE_TOOLS_DEFINED_H

#include <Eigen/Core>

namespace ACE{

extern Eigen::VectorXcd H_Matrix_to_L_Vector(const Eigen::MatrixXcd &M);

extern Eigen::MatrixXcd L_Vector_to_H_Matrix(const Eigen::VectorXcd &v);

/* Outer product for operators in Liouville space: 
(alpha1, alpha2) otimes (beta1, beta2) -> (gamma1, gamma2)
((nu1,mu1),(nu2,mu2)) otimes ((xi1,chi1),(xi2,chi2)) 
-> ((nu1 *dim(xi1) + xi1), mu1*dim(chi1) + chi1) , ((nu2* .. 
*/
extern Eigen::MatrixXcd L_Op_otimes(const Eigen::MatrixXcd & A, 
                                   const Eigen::MatrixXcd & B);

/* For propagator:  Hilbert_total = Hilbert_S otimes Hilbert_E
                    Liouville_total = Hilbert_total otimes Hilbert_total
   disentangle as:  Liouville_total = Liouville_S otimes Liouville_E
*/ 
extern Eigen::MatrixXcd Disentangle_Propagator(const Eigen::MatrixXcd &M, int Ns);

//Expand Hilbert space of matrix: 
//HS otimes HE -> HS otimes 1_by_dim otimes HE
//with d1=dim(HS)
extern Eigen::MatrixXcd ExpandMatrix(const Eigen::MatrixXcd &A, int d1, int by_dim);

// Same for operators in Liouville space. 
// d1 and by_dim are still Hilbert space dimensions
extern Eigen::MatrixXcd ExpandMatrixLiouville(const Eigen::MatrixXcd &A, int d1, int by_dim);

extern bool print_diff_from_ortho(const Eigen::MatrixXcd &A, double threshold=1e-12, const std::string &str="");

extern bool print_diff_from_ortho(const Eigen::MatrixXd &A, double threshold=1e-12, const std::string &str="");

extern double max_diff_from_ortho(const Eigen::MatrixXcd &A);

extern double max_diff_from_ortho(const Eigen::MatrixXd &A);

extern void print_nonzero(const Eigen::MatrixXcd &A, double epsilon=1e-20);

}//namespace
#endif
