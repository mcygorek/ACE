#ifndef INFLUENCE_FUNCTIONAL_DEFINED_H
#define INFLUENCE_FUNCTIONAL_DEFINED_H

#include <Eigen/Core>
#include "SpectralDensity.h"
#include "Tensor.h"
#include "DiagBB.h"
#include <fstream>

class InfluenceFunctional{
public:
  int n_max;
  double dt;

  DiagBB diagBB;
  std::vector<Tensor_Dense> ten;
  Coupling_Groups groups;

  inline int get_grp(int i)const{return groups.grp[i];}  
  inline int get_Ngrps()const{return groups.Ngrps;}

  virtual int get_dim()const{return diagBB.get_dim();}
//  virtual int get_dim()const{return couplings.size();}

  int get_n_max()const{return n_max;}
  double get_dt()const{return dt;}

  Tensor_Dense & operator[](int i){ return ten[i];}
  const Tensor_Dense & operator[](int i)const{ return ten[i];}

  void calculate(){
    int NL=get_dim()*get_dim();

#ifdef DEBUG
std::cout<<"InfluenceFunctional::calculated called with n_max="<<n_max<<" and NL="<<NL<<std::endl;
{int bs=1; for(int i=0; i<n_max+1; i++)bs*=NL;
std::cout<<"InfluenceFunctional: Total blocksize: "<<bs<<std::endl;}
#endif

    
    ten.resize(n_max+1);
    std::vector<Eigen::MatrixXcd> eS(n_max+1);
    for(int i=0; i<n_max+1; i++){
      ten[i].resize(i+1, NL);
      eS[i]=diagBB.calculate_expS(i, get_dt());


      if(i==0){
        for(int j=0; j<NL; j++){
          ten[0][j]=eS[0](j,j);
        }
      }else{
        int bs=ten[i-1].get_total_size()/NL;
        for(int b=0; b<bs; b++){
          for(size_t k=0; k<NL; k++){
            for(size_t l=0; l<NL; l++){
              ten[i][(k*bs+b)*NL+l]=ten[i-1][k*bs+b]*eS[i](k,l);
            }
          }
        }
      }
    }
  }


  void print(const std::string &fname)const{
    std::ofstream ofs(fname.c_str());
    ofs<<"n_max: "<<n_max<<std::endl;
    for(size_t n=0; n<ten.size(); n++){
      for(Tensor_Index ind(ten[n]); !ind.done(ten[n]); ind.increment(ten[n])){

        ofs<<"n: "<<n<<" index: "<<ind<<" value: "<<ten[n](ind)<<std::endl;
      }
    }
  }
  InfluenceFunctional(int n_max_, double dt_, 
                      const Eigen::MatrixXcd &couplings_,
                      RealFunctionPtr SD_, double temperature_, 
                      bool noSubPS=false)
   : n_max(n_max_), dt(dt_), groups(couplings_),
     diagBB(Coupling_Groups(couplings_), couplings_, SD_, temperature_, noSubPS) {

    calculate();
  }
  ~InfluenceFunctional(){
  }
};

#endif
