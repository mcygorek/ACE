#ifndef ACE_DYNAMICALMAP_DEFINED_H
#define ACE_DYNAMICALMAP_DEFINED_H

#include "Simulation_PT.hpp"
#include "LindbladMasterEquation.hpp"

namespace ACE{

/**
Calculates and stores time-dependent dynamical maps, i.e. the E in 
\bar{\rho}(t)= E_{t,t_0} \bar{\rho}(t_0)
where \bar{\rho}(t) is the reduced density matrix.
*/


class DynamicalMap{
public:
double ta, dt;
std::vector<Eigen::MatrixXcd> E;
std::vector<Eigen::MatrixXcd> Pade;

const Eigen::MatrixXcd &get(int i)const;
static Eigen::MatrixXcd invert(const Eigen::MatrixXcd &M, double regularize=0., std::ostream *os=NULL);
Eigen::MatrixXcd get_dE(int i, double regularize=0., std::ostream *os=NULL)const;
Eigen::MatrixXcd get_dE_ref(int i, double eps, double tref)const;


void calculate(Propagator &prop, ProcessTensorForwardList &PT, Simulation_PT &sim, const TimeGrid &tgrid, int verbosity=0);

void calculate(Parameters &param);
void calculate_TEMPO(Parameters &param);
void calculate_stepwise(Parameters &param);

static void make_physical_dE(Eigen::MatrixXcd & dE);
void make_physical();

void propagate(const TimeGrid &tgrid, const Eigen::MatrixXcd &intial_state,  OutputPrinter &printer)const;
void propagate_dE(const TimeGrid &tgrid, const Eigen::MatrixXcd &initial_state,  OutputPrinter &printer)const;

Eigen::MatrixXcd get_Liouvillian(int l, double regularize=0.)const;
LindbladMasterEquation get_LindbladMasterEquation(int l)const;
LindbladMasterEquation get_LindbladMasterEquation_Hall(int l)const;
void print_Lindblad_Hall(const std::string & prefix)const;
void print_Lindblad(const std::string & prefix)const;


void print_normdiff(std::ostream &os)const;
inline void print_normdiff(const std::string &fname)const{
  std::ofstream ofs(fname.c_str());
  print_normdiff(ofs);
}
void print_TT_normdiff(std::ostream &os)const;
inline void print_TT_normdiff(const std::string &fname)const{
  std::ofstream ofs(fname.c_str());
  print_TT_normdiff(ofs);
}
void print_TT_norm(std::ostream &os)const;
inline void print_TT_norm(const std::string &fname)const{
  std::ofstream ofs(fname.c_str());
  print_TT_norm(ofs);
}

void print_eigenvalues(std::ostream &os, double regularize=0.)const;
inline void print_eigenvalues(const std::string &fname, double regularize=0.)const{
  std::ofstream ofs(fname.c_str());
  print_eigenvalues(ofs,regularize);
}
void print_dE_eigenvalues(std::ostream &os)const;
inline void print_dE_eigenvalues(const std::string &fname)const{
  std::ofstream ofs(fname.c_str());
  print_dE_eigenvalues(ofs);
}
void print_E_eigenvalues(std::ostream &os)const;
inline void print_E_eigenvalues(const std::string &fname)const{
  std::ofstream ofs(fname.c_str());
  print_E_eigenvalues(ofs);
}
void print_eigenvalues_ref(std::ostream &os, double eps, double tref)const;
inline void print_eigenvalues_ref(const std::string &fname, double eps, double tref)const{
  std::ofstream ofs(fname.c_str());
  print_eigenvalues_ref(ofs,eps,tref);
}
void print_dE_eigenvalues_ref(std::ostream &os, double eps, double tref)const;
inline void print_dE_eigenvalues_ref(const std::string &fname, double eps, double tref)const{
  std::ofstream ofs(fname.c_str());
  print_dE_eigenvalues_ref(ofs,eps,tref);
}

void print_fullE_eigenvalues(std::ostream &os)const;
inline void print_fullE_eigenvalues(const std::string &fname)const{
  std::ofstream ofs(fname.c_str());
  print_fullE_eigenvalues(ofs);
}
void print_dE_singularvalues(std::ostream &os)const;
inline void print_dE_singularvalues(const std::string &fname)const{
  std::ofstream ofs(fname.c_str());
  print_dE_singularvalues(ofs);
}
void print_singularvalues(std::ostream &os)const;
inline void print_singularvalues(const std::string &fname)const{
  std::ofstream ofs(fname.c_str());
  print_singularvalues(ofs);
}


Eigen::MatrixXcd get_average_dE(double tstart, double tend)const;
void Richardson_combine(const DynamicalMap &other);


void print_average_eigenvalues(std::ostream &os, double tstart, double tend)const;
inline void print_average_eigenvalues(const std::string &fname, double tstart, double tend)const{
  std::ofstream ofs(fname.c_str());
  print_average_eigenvalues(ofs, tstart, tend);
}

void analyze_physicality(std::ostream &os)const;
inline void analyze_physicality(const std::string &fname)const{
  std::ofstream ofs(fname.c_str());
  analyze_physicality(ofs);
}

std::vector<Eigen::MatrixXcd> get_TT(int length=0)const;
void extrapolate_TT(double te2, int length=0);

void calculate_Pade(const std::vector<int> &steps);
void calculate_Pade2(const std::vector<int> &steps);
void Pade_extrapolate(int add_steps);

void read(const std::string &fname);
void write(const std::string &fname)const;


DynamicalMap(){}
DynamicalMap(const std::string &fname){ read(fname); }
~DynamicalMap(){}

};
}//namespace
#endif
