#ifndef INFLUENCE_FUNCTIONAL_MPS_DEFINED_H
#define INFLUENCE_FUNCTIONAL_MPS_DEFINED_H

#include <Eigen/Core>
#include "SpectralDensity.h"
#include "MPS.h"
#include "DiagBB.h"
#include <fstream>

class InfluenceFunctional_Vector{
public:
  double dt;

  DiagBB diagBB;
  Coupling_Groups groups;
  Coupling_Groups_Liouville lgroups;
 
  std::vector<Eigen::MatrixXcd> b;

  inline int get_grp(int i)const{return groups.grp[i];}  
  inline int get_Ngrps()const{return groups.Ngrps;}

  inline int get_grp2(int i)const{return lgroups.grp[i];}  
  inline int get_Ngrps2()const{return lgroups.Ngrps;}

//  virtual int get_dim()const{return diagBB.get_dim();}
  virtual int get_dim()const{return groups.sys_dim();}

  int get_n_max()const{return b.size()-1;}
  double get_dt()const{return dt;}


  void calculate(int n_max){
    int NL=get_dim()*get_dim();
    if(n_max<0){
      std::cerr<<"n_max must not be negative!"<<std::endl;
      exit(1);
    }
    
    b.resize(n_max+1);
    for(int n=0; n<n_max+1; n++){
//      b[n]=diagBB.calculate_expS(n_max-n, get_dt());
      b[n]=diagBB.calculate_expS(n, get_dt());
    }   
  }


  void print(const std::string &fname)const{
    std::cerr<<"InfluenceFunctional_Vector: print not implemented yet!"<<std::endl;
    exit(1);
  }
  InfluenceFunctional_Vector(int n_max_, double dt_, 
                      const Eigen::MatrixXcd &couplings_,
                      RealFunctionPtr SD_, double temperature_)
   : dt(dt_), groups(couplings_), lgroups(couplings_),
     diagBB(Coupling_Groups(couplings_), couplings_, SD_, temperature_) {


    calculate(n_max_);
  }
  ~InfluenceFunctional_Vector(){
  }
};

#endif
