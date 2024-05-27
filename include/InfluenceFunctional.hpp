#ifndef INFLUENCE_FUNCTIONAL_DEFINED_H
#define INFLUENCE_FUNCTIONAL_DEFINED_H

#include <vector>
#include <Eigen/Core>
#include "SpectralDensity.hpp"
#include "Tensor_Dense.hpp"
#include "DiagBB.hpp"
#include <fstream>

namespace ACE{

class InfluenceFunctional{
public:
  int n_max;
  double dt;

  std::vector<Tensor_Dense> ten;
  Coupling_Groups groups;
  DiagBB diagBB;

  inline int get_grp(int i)const{return groups.grp[i];}  
  inline int get_Ngrps()const{return groups.Ngrps;}

  inline virtual int get_dim()const{return diagBB.get_dim();}

  inline int get_n_max()const{return n_max;}
  inline double get_dt()const{return dt;}

  inline Tensor_Dense & operator[](int i){ return ten[i];}
  inline const Tensor_Dense & operator[](int i)const{ return ten[i];}

  void calculate();

  void print(const std::string &fname)const;
  
  void setup(int n_max_, double dt_, DiagBB &diagBB_);
  
  inline InfluenceFunctional(int n_max_, double dt_, 
                      const Eigen::MatrixXcd &couplings_,
                      RealFunctionPtr SD_, double temperature_, 
                      bool noSubPS=false, 
                      double omega_min=0., double omega_max=50., 
                      double E_shift_init=0.)
   : n_max(n_max_), dt(dt_), groups(couplings_),
     diagBB(Coupling_Groups(couplings_), couplings_, SD_, temperature_, noSubPS,omega_min, omega_max, E_shift_init) {

    calculate();
  }
  inline InfluenceFunctional(int n_max_, double dt_, DiagBB &diagBB_)
   : n_max(n_max_), dt(dt_), diagBB(diagBB_) {
    calculate();
  }
  inline InfluenceFunctional(){
  }
  inline ~InfluenceFunctional(){
  }
};

}//namespace
#endif
