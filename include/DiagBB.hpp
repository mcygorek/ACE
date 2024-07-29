#ifndef DIAG_BB_DEFINED_H
#define DIAG_BB_DEFINED_H

#include <vector>
#include "Function.hpp"
#include "Coupling_Groups.hpp"
#include "HilbertSpaceRotation.hpp"

namespace ACE{
class Parameters;

class DiagBB{
public:
  RealFunctionPtr J;
  double temperature;
  double E_shift_init;
  double omega_min;
  double omega_max;
  double damping;
  double damping_Gaussian;
  size_t Ndiscr;
  bool noSubPS;
  bool separate_freq;
  bool high_T_limit;
  HilbertSpaceRotation hs_rot;
  std::vector<double> couplings;
  Coupling_Groups groups;

  std::vector<std::complex<double> > K_precalc;
  
  inline bool is_set_up()const{ return groups.is_set_up(); }
  virtual int get_dim()const{ return couplings.size(); }
  virtual int sys_dim()const;

  double get_beta()const; 
  static double get_coth(double beta, double E_shift, double w);

#include "DiagBB_K_integrands.hpp"


  //obtain K_XX for diagonally coupled bath
  std::complex<double> calculate_K(int n, double dt);
  std::complex<double> calculate_K_explicit(int n, double dt);

  Eigen::MatrixXcd calculate_expS(int n, double dt);

  void print_K(const std::string &fname, int n_max, double dt);

  void read_K_int(const std::string &fname, int n_max, double dt);
 
  void precalc_FFT(int n_max, double dt);

    // Find memory time: relative to maximal (max) value around K(t=0).
    // K can be oscillatory. So, after finding first |K| < threshold * max, go at least 2 twice as long to see if values larger than that are found
  int estimate_memory_length(int n_max, double dt, double threshold, bool verbose);
  
  void setup_groups_and_couplings(const Eigen::MatrixXcd &Op);

  void setup(const Coupling_Groups & grp, 
                const Eigen::MatrixXcd &couplings_, 
                RealFunctionPtr J_, double temperature_=0., 
                bool noSubPS_=false,
                double omega_min_ = 0., double omega_max_ = 50., 
                double E_shift_init_ =0., bool high_T=false );
  
  static void complain_if_not_Hermitian(const Eigen::MatrixXcd & sysop);
  
  static bool is_offdiagonal(const Eigen::MatrixXcd & sysop);
  
  void setup(Parameters &param, const std::string &prefix);

  inline DiagBB(const Coupling_Groups & grp, 
                const Eigen::MatrixXcd &couplings_, 
                RealFunctionPtr J_, double temperature_=0., 
                bool noSubPS_=false,
                double omega_min_ = 0., double omega_max_ = 50., 
                double E_shift_init_ =0., bool high_T=false ) : J(J_){

    setup(grp, couplings_, J_, temperature_, noSubPS_, omega_min_, omega_max_, E_shift_init_, high_T);
  }
  inline DiagBB(Parameters &param, const std::string &prefix){
    setup(param, prefix);
  }
  DiagBB() : J (RealFunctionPtr_Zero) {}
};

}//namespace
#endif
