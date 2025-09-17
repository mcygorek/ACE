#ifndef INFLUENCE_FUNCTIONAL_VECTOR_DEFINED_H
#define INFLUENCE_FUNCTIONAL_VECTOR_DEFINED_H

#include <vector>
#include <Eigen/Core>
#include "SpectralDensity.hpp"
#include "MPS.hpp"
#include "DiagBB.hpp"

namespace ACE{

class InfluenceFunctional_Vector{
public:
  double dt;

  Coupling_Groups groups;
  Coupling_Groups_Liouville lgroups;
  DiagBB diagBB;
 
  std::vector<Eigen::MatrixXcd> b;

  inline int get_grp(int i)const{return groups.grp[i];}  
  inline int get_Ngrps()const{return groups.Ngrps;}

  inline int get_grp2(int i)const{return lgroups.grp[i];}  
  inline int get_Ngrps2()const{return lgroups.Ngrps;}

//  virtual int get_dim()const{return diagBB.get_dim();}
  inline virtual int get_dim()const{return groups.sys_dim();}

  inline int get_n_max()const{return b.size();}
  inline double get_dt()const{return dt;}


  void calculate(int n_max);

  void print(const std::string &fname)const;
  
  inline InfluenceFunctional_Vector(int n_max_, double dt_, 
                      const Eigen::MatrixXcd &couplings_,
                      RealFunctionPtr SD_, double temperature_,
                      bool noSubPS=false,
                      double omega_min=0., double omega_max=50.,
                      double E_shift_init=0.)
   : dt(dt_), groups(couplings_), lgroups(couplings_),
     diagBB(Coupling_Groups(couplings_), couplings_, SD_, temperature_, noSubPS, omega_min, omega_max, E_shift_init) {


    calculate(n_max_);
  }
  inline void setup(int n_max_, double dt_, const DiagBB &diagBB_){
    dt=dt_;
    groups=diagBB_.groups;
    lgroups=diagBB_.groups;
    diagBB=diagBB_;

    calculate(n_max_);
  }
  inline InfluenceFunctional_Vector(int n_max_, double dt_, const DiagBB &diagBB_)
   : dt(dt_), groups(diagBB_.groups), lgroups(diagBB_.groups),
     diagBB(diagBB_) { 

    calculate(n_max_);
  }
  inline InfluenceFunctional_Vector(){
  }
  inline ~InfluenceFunctional_Vector(){
  }
};

}//namespace
#endif
