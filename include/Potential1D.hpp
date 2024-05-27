#ifndef ACE_POTENTIAL1D_DEFINED_H
#define ACE_POTENTIAL1D_DEFINED_H

#include <vector>
#include <Eigen/Dense>

/**
Solution of Schrödinger equation of a one-dimensional problem. 
-> extraction of effective model based on the lowest M eigenstates 

*/
namespace ACE{
class Parameters;

class Potential1D{
public:

  int M; //Nr eigenstates extracted

  std::vector<double> vx; //unitless potential discretized on regular grid
  double xa, dx; //defines grid (together with vx.size())

  double E_scale, x_scale; //scaling of units
  double ddx2_scale; //prefactor in front of d2^/dx^2 in ODE
  double x_coupling; // factor to link (a^\dagger + a) = x_coupling * <n|x|m>  


  Eigen::VectorXd E; //result: eigenvalues
  Eigen::MatrixXd X; //result: matrix elements: position operator
  Eigen::MatrixXd X2; //result: matrix elements: position operator spread


  void calculate(const std::string &printEV="");

  std::string add_name(const std::string &name, const std::string &str);

  void setup(Parameters &param, const std::string &name="");

  inline Potential1D(Parameters &param, const std::string &name=""){
    setup(param,name);
  }
  Potential1D();
};

}//namespace
#endif
