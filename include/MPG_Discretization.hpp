#ifndef MPG_DISCRETIZATION_DEFINED_H
#define MPG_DISCRETIZATION_DEFINED_H

#include <vector>
#include <iostream>
#include "EnergyRange.hpp"
#include "SpectralDensity_Selector.hpp"

namespace ACE{

class MPG_Discretization{
public:
  int N;
  std::vector<double> E;
  std::vector<double> dE;

  virtual void check_bounds(const std::string fct, int k)const;

  double get_E(int k)const;
  double get_dE(int k)const;
  double get_omega(int k)const;
  double get_domega(int k)const;

  std::string add_name(const std::string &mpgname, const std::string &str)const;

  virtual void setup(Parameters &param, const std::string &mpgname);

  inline MPG_Discretization(Parameters &param, const std::string &mpgname){
    setup(param, mpgname);
  }
  inline MPG_Discretization(){}
};


//MPG_Discretization for a one-parameter interaction
class MPG_Discretization_E_g: public MPG_Discretization{
public:
  std::vector<double> g;

  void zero_pad(int N_new);

  virtual void check_bounds(const std::string fct, int k)const;
  
  double get_g(int k)const;

  static double rate_from_g(double g, double DOS);
  static double g_from_rate(double rate, double DOS);

  static bool compare_abs(const std::pair<double,int> &p1,
                          const std::pair<double,int> &p2);
  static bool compare_abs_less(const std::pair<double,int> &p1,
                          const std::pair<double,int> &p2);
  static bool compare_less(const std::pair<double,int> &p1,
                          const std::pair<double,int> &p2);
  static bool compare_greater(const std::pair<double,int> &p1,
                          const std::pair<double,int> &p2);

  void print(std::ostream &os=std::cout)const;
  void print(const std::string &fname)const;
  
  virtual void setup(Parameters &param, const std::string &mpgname);

  inline MPG_Discretization_E_g(Parameters &param, const std::string &mpgname=""){
    setup(param, mpgname);
  }
  inline MPG_Discretization_E_g(){}
};

}//namespace
#endif
