#include "ADM.hpp"
#include "RankCompressor.hpp"
#include "Propagator.hpp"
#include "InfluenceFunctional.hpp"
#include "Tensor_Dense.hpp"
#include <Eigen/Dense>
#include "LiouvilleTools.hpp"

namespace ACE{

  void AugmentedDensityMatrix::check_dimensions()const{
    if(get_n_max()<1){
      std::cerr<<"AugmentedDensityMatrix::check_dimensions: get_n_max()<1!"<<std::endl; 
      exit(1);
    }
    int N=sqrt(get_NL());
    if(N*N != get_NL()){
      std::cerr<<"AugmentedDensityMatrix::check_dimensions: N*N!=NL!"<<std::endl; 
      exit(1);
    }
  }
  
  void AugmentedDensityMatrix::print_status(std::ostream &os)const{
  }

  void AugmentedDensityMatrix::update_rho(int step){
    //Note: initialization (step 0): n_mem_eff=1
    //      step 1: after special contraction n_mem_eff=1
    //      step 2: after expansion n_mem_eff=2 ...
    int n_mem_eff=step;
    if(step<1)n_mem_eff=1;
    if(n_mem_eff>get_n_max())n_mem_eff=get_n_max();
    int NL=get_NL();
    int N=sqrt(NL);

    int Ngrps2=get_Ngrps2();

    int back_blocksize=1;
    for(int l=1; l<n_mem_eff; l++)back_blocksize*=Ngrps2;

    rho=Eigen::MatrixXcd::Zero(N,N);
    for(int i=0; i<N; i++){
      for(int j=0; j<N; j++){
        for(int b=0; b<back_blocksize; b++){
          rho(i,j)+=ten[(i*N+j)*back_blocksize+b];
        }
      }
    }
  }

  void AugmentedDensityMatrix::propagate(
                 Propagator &prop, const InfluenceFunctional &IF, 
                 double t, double dt, int step, RankCompressor *compressor,
                 bool use_symmetric_Trotter){
    
//It's not completely straightforward to define a QUAPI algorithm with symmetric Trotter splitting. This is because we have to connect to past points between time steps. We solve this by applying only a half-step M only in the first time step and use full steps for all other steps. The missing half-steps are added after update_rho(..). This implies that one should never call update_rho(..) from any other function than this one.


//    if(use_symmetric_Trotter){
//      std::cout<<"--------------------------------------------------------------------------"<<std::endl;
//      std::cout<<"WARNING: using QUAPI with use_symmetric_Trotter=true might be problematic!"<<std::endl;
//      std::cout<<"--------------------------------------------------------------------------"<<std::endl;
//    }

    check_dimensions();
    int n_mem=get_n_max()+1;
    int NL=get_NL();
    int N=sqrt(NL);

    int Ngrps=IF.get_Ngrps();
    int Ngrps2=get_Ngrps2();
    if(Ngrps*Ngrps!=Ngrps2){
      std::cerr<<"ADM::propagate: Ngrps*Ngrps!=Ngrps2"<<std::endl;
      std::cerr<<"IF.get_Ngrps(): "<<IF.get_Ngrps()<<" Ngrps2: "<<get_Ngrps2()<<std::endl;
      exit(1);
    }

    std::vector<int> grp2(NL);
    for(int i=0; i<N; i++){
      for(int j=0; j<N; j++){
        grp2[i*N+j]=IF.get_grp(i)*Ngrps+IF.get_grp(j);
      }
    }

    Tensor_Dense ten2(ten.get_dims());
    ten2.fill(0);

    if(use_symmetric_Trotter){
      if(step == 1){ 
        prop.update(t,dt/2.);
      }else{
        prop.update(t-dt/2.,dt);
      }
    }else{
      prop.update(t,dt);
    }

    if(step <1){
      std::cerr<<"AugmentedDensityMatrix::propagate: step <1!"<<std::endl;
      exit(1);
    }else if(step == 1 || n_mem < 2){//first step: contract last index
      for(int i=0; i<NL; i++){
        for(int j=0; j<NL; j++){
          ten2[i] += IF.ten[0][grp2[i]]* prop.M(i,j) * ten[j];
        }
      }
    }else if(step<n_mem){ //step with expansion
      int BS=1;
      for(int i=0; i<step-2; i++)BS*=Ngrps2;

      for(int i=0; i<NL; i++){
        int gi=grp2[i];
        for(int j=0; j<NL; j++){ 
          int gj=grp2[j];
          for(int b=0; b<BS; b++){
            ten2[(i*Ngrps2+gj)*BS+b] +=
              IF.ten[step-1][(gi*Ngrps2+gj)*BS+b] * 
              prop.M(i,j) * ten[j*BS+b];
          }
        }
      }
    }else if(step==n_mem && n_mem==2){ //M and truncation over same index
      for(int i=0; i<NL; i++){
        int gi=grp2[i];
        for(int j=0; j<NL; j++){ 
          int gj=grp2[j];
          ten2[i] += IF.ten[step-1][gi*Ngrps2+gj] * prop.M(i,j) * ten[j];
        }
      }
    }else{ //propagate and contract last index   
      int BS=1;
      for(int i=0; i<n_mem-3; i++)BS*=Ngrps2;

      for(int i=0; i<NL; i++){
        int gi=grp2[i];
        for(int j=0; j<NL; j++){ 
          int gj=grp2[j];
          for(int b=0; b<BS; b++){
            for(int l=0; l<Ngrps2; l++){ 
              ten2[(i*Ngrps2+gj)*BS+b] +=
                IF.ten[n_mem-1][((gi*Ngrps2+gj)*BS+b)*Ngrps2+l] * 
                prop.M(i,j) * ten[(j*BS+b)*Ngrps2+l];
            }
          }
        }
      }

    }

    ten.swap(ten2);

    update_rho(step);
    if(use_symmetric_Trotter){
      prop.update(t+dt/2.,dt/2.);
      rho=L_Vector_to_H_Matrix( prop.M * H_Matrix_to_L_Vector(rho) );
    }
  }
  
  
  AugmentedDensityMatrix::AugmentedDensityMatrix(
            int n_max, int Ngrps, const Eigen::MatrixXcd &rho_)  : rho(rho_) {
    int N=rho_.rows();
    if(rho.rows()!=rho.cols() ){
      std::cerr<<"ADM: constructor: rho has wrong dimensions!"<<std::endl;
      exit(1);
    }

    std::vector<int> tensordims(n_max, Ngrps*Ngrps);
    if(n_max>0)tensordims[0]=N*N;
    ten.resize(tensordims);

    if(n_max>0){
      ten.fill(0.);
      for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
          ten[i*N+j]=rho(i,j);
        }
      }
    }
  }

}
