#include "InfluenceFunctional.hpp"
#include <Eigen/Core>
#include "SpectralDensity.hpp"
#include "Tensor_Dense.hpp"
#include "DiagBB.hpp"
#include <fstream>

namespace ACE{

  void InfluenceFunctional::calculate(){
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
          for(int k=0; k<NL; k++){
            for(int l=0; l<NL; l++){
              ten[i][(k*bs+b)*NL+l]=ten[i-1][k*bs+b]*eS[i](k,l);
            }
          }
        }
      }
    }
  }


  void InfluenceFunctional::print(const std::string &fname)const{
    std::ofstream ofs(fname.c_str());
    ofs<<"n_max: "<<n_max<<std::endl;
    for(size_t n=0; n<ten.size(); n++){
      for(Tensor_Index ind(ten[n]); !ind.done(ten[n]); ind.increment(ten[n])){

        ofs<<"n: "<<n<<" index: "<<ind<<" value: "<<ten[n](ind)<<std::endl;
      }
    }
  }

  void InfluenceFunctional::setup(int n_max_, double dt_, DiagBB &diagBB_){
    n_max=n_max_;
    dt=dt_;  
    diagBB=diagBB_;
    calculate();
  }

}//namespace
