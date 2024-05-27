#ifndef MODE_PROPAGATOR_GENERATOR_RANDOMSPIN_DEFINED_H
#define MODE_PROPAGATOR_GENERATOR_RANDOMSPIN_DEFINED_H

#include "ModePropagatorGenerator.hpp"
#include <vector>
#include <Eigen/Dense>

namespace ACE{

class ModePropagatorGenerator_RandomSpin: public ModePropagatorGenerator{
public:

  std::vector<double> J;  
  std::vector<Eigen::Vector3d> initialDirs;
  std::vector<Eigen::Vector3d> B_eff;

  inline double get_J(int k)const{ return J[k]; }

  inline virtual std::string name()const{return std::string("RandomSpin");}

  inline virtual int get_N()const{ return 2; }

  static bool compare_abs_smaller(const double &p1, const double &p2);
  static bool compare_abs_larger(const double &p1, const double &p2);
  static bool compare_dir_x(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2);
  static bool compare_dir_z(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2);

  virtual std::vector<Eigen::MatrixXcd> get_env_ops(int k)const;

  void setup(int Nmod, double J_max, double J_min=0, size_t seed=1);

  virtual void setup(Parameters &param);

  virtual ModePropagatorPtr getModePropagator(int k)const;
  
  virtual Eigen::MatrixXcd get_bath_init(int k)const;

  inline ModePropagatorGenerator_RandomSpin(int Nmod, double Jmax, double Jmin=0){
    setup(Nmod, Jmax, Jmin);
  }
  inline ModePropagatorGenerator_RandomSpin(Parameters &param){
    setup(param);
  }
};

}//namespace
#endif
