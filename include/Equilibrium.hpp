#ifndef EQUILIBRIUM_DEFINED_H
#define EQUILIBRIUM_DEFINED_H

#include "Eigen_fwd.hpp"

namespace ACE{

Eigen::MatrixXcd Nlevel_Equilibrium(const Eigen::VectorXcd &E, double T);

Eigen::MatrixXcd Boson_Equilibrium(int n_max, double dE, double T);

// Equilibrium w.r.t. to Hamiltonian matrix
Eigen::MatrixXcd Boson_Equilibrium(Eigen::MatrixXcd H, double T, double E_shift=0.);


//calculates the factor x=dE/(kB*T), so that equilibrium with Boltzmann factor e^{-x n} yields the average boson number 'nav'
double Boson_Equilibrium_invert(int nmax, double nav);

}//namespace
#endif
