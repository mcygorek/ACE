#pragma once
#ifndef ACE_LINDBLAD_MASTER_EQUATION_DEFINED_H
#define ACE_LINDBLAD_MASTER_EQUATION_DEFINED_H

#include <iostream>
#include <Eigen/Eigenvalues>
#include <vector>
#include "TimeGrid.hpp"

namespace ACE{

struct LindbladMasterEquation{
  Eigen::MatrixXcd H;                           // Hamiltonian
  std::vector<std::pair<double,Eigen::MatrixXcd> > L; // (rate, operator)

  int get_dim()const;

  Eigen::MatrixXcd construct_Liouvillian()const;
  
  void set_from_Liouvillian_Hall(const Eigen::MatrixXcd &Liou, double threshold, int verbosity);
  void set_from_Liouvillian(const Eigen::MatrixXcd &Liou, double threshold, int verbosity);

  void print(std::ostream &os=std::cout, double epsilon=0)const;
  void print_param(const std::string &fname, const TimeGrid &tgrid, bool print_Markov, double epsilon)const;

  LindbladMasterEquation(){}
  LindbladMasterEquation(const Eigen::MatrixXcd &Liou, double threshold, int verbosity){
    set_from_Liouvillian(Liou, threshold, verbosity);
  }
};
}//namespace


#endif
