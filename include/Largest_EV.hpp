#ifndef ACE_LARGEST_EV_DEFINED_H
#define ACE_LARGEST_EV_DEFINED_H

#include "Eigen_fwd.hpp"
#include <functional>
#include <complex>

namespace ACE{

std::complex<double> Largest_EV_Arnoldi(Eigen::VectorXcd &initial, int m, 
      double epsilon,
      std::function<Eigen::VectorXcd(const Eigen::VectorXcd &v)> Afunc, 
      int verbosity, bool reortho=false);

std::complex<double> Largest_EV_Arnoldi_BLAS(Eigen::VectorXcd &initial, int m, 
      double epsilon,
      std::function<Eigen::VectorXcd(const Eigen::VectorXcd &v)> Afunc, 
      int verbosity);

std::complex<double> Largest_EV_KrylovSchur(Eigen::VectorXcd &vec, 
      int itertotal, int m, int k, double epsilon,
      std::function<Eigen::VectorXcd(const Eigen::VectorXcd &v)> Afunc, 
      int verbosity);

std::complex<double> Largest_EV_Arnoldi_restart(Eigen::VectorXcd &vec, 
      int maxiter, double epsilon, 
      std::function<Eigen::VectorXcd(const Eigen::VectorXcd &v)> Afunc,
      int verbosity, bool reortho=false, int m=10);

}//namespace
#endif
