#ifndef INFLUENCE_FUNCTIONAL_MPS_DEFINED_H
#define INFLUENCE_FUNCTIONAL_MPS_DEFINED_H

#include <Eigen/Core>
#include "SpectralDensity.h"
#include "MPS.h"
#include "DiagBB.h"
#include <fstream>

class InfluenceFunctional_MPS: public MPS{
public:
  double dt;

  DiagBB diagBB;
  Coupling_Groups groups;

  inline int get_grp(int i)const{return groups.grp[i];}  
  inline int get_Ngrps()const{return groups.Ngrps;}

  virtual int get_dim()const{return diagBB.get_dim();}
//  virtual int get_dim()const{return couplings.size();}

  int get_n_max()const{return MPS::get_rank()-1;}
  double get_dt()const{return dt;}


  void calculate(int n_max){
    int NL=get_dim()*get_dim();
    if(n_max<0){
      std::cerr<<"n_max must not be negative!"<<std::endl;
      exit(1);
    }
    
    a.resize(n_max+1);
    if(n_max<1){
      Eigen::MatrixXcd eS=diagBB.calculate_expS(0, get_dt());
      a[0].resize(NL,1,1);
      for(int i=0; i<NL; i++){
        a[0](i,0,0)=eS(i,i);
      }      
      return;
    }


    //first
    a[0].resize(NL,1,NL);
    Eigen::MatrixXcd eS=diagBB.calculate_expS(n_max, get_dt());
    for(int i=0; i<NL; i++){
      for(int j=0; j<NL; j++){
        a[0](i,0,j)=eS(i,j);
      }
    }
    //last
    a[n_max].resize(NL,NL,1);
    eS=diagBB.calculate_expS(0, get_dt());
    for(int i=0; i<NL; i++){
      for(int j=0; j<NL; j++){
        a[0](i,j,0)=0.;
        if(j==i)a[0](i,j,0)=eS(i,i);
      }
    }
  
    //rest

    for(int n=1; n<n_max; n++){
      eS=diagBB.calculate_expS(n_max-n, get_dt());
      a[n].resize(NL,NL,NL);
      a[n].fill(0.);
      for(int i=0; i<NL; i++){
        for(int j=0; j<NL; j++){
          a[n](i,j,j)=eS(i,j);
        }
      }
    }   
  }


  void print(const std::string &fname)const{
    std::ofstream ofs(fname.c_str());
    int n_max=get_n_max();
    ofs<<"n_max: "<<n_max<<std::endl;
    for(Tensor_Index ind(*this); !ind.done(*this); ind.increment(*this)){
      ofs<<"n: "<<n_max<<" index: "<<ind<<" value: "<<MPS::operator()(ind)<<std::endl;
    }
  }
  InfluenceFunctional_MPS(int n_max_, double dt_, 
                      const Eigen::MatrixXcd &couplings_,
                      RealFunctionPtr SD_, double temperature_,
                      bool noSubPS=false)
   : dt(dt_), groups(couplings_),
     diagBB(Coupling_Groups(couplings_), couplings_, SD_, temperature_,noSubPS) {


    calculate(n_max_);
  }
  ~InfluenceFunctional_MPS(){
  }
};

#endif
