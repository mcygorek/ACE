#include "ADM.hpp"
#include "RankCompressor.hpp"
#include "Propagator.hpp"
#include "InfluenceFunctional.hpp"
#include "Tensor_Dense.hpp"
#include <Eigen/Dense>

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

    check_dimensions();
    int n_mem=get_n_max();
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

//std::cout<<"NL: "<<NL<<" n_mem: "<<n_mem<<" Ngrps: "<<Ngrps<<" ";
//std::cout<<"ADM dims: "; ten.print_dims(std::cout); std::cout<<std::endl;

    Tensor_Dense ten2(ten.get_dims());
    ten2.fill(0);

    Eigen::MatrixXcd M2=Eigen::MatrixXcd::Identity(NL,NL);
    if(use_symmetric_Trotter){
      prop.update(t+dt/2.,dt/2.);
      M2=prop.M;
      prop.update(t,dt/2.);
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
      ten.swap(ten2);
      ten2.fill(0.);
      if(use_symmetric_Trotter){
        for(int i=0; i<NL; i++){
          for(int j=0; j<NL; j++){
            ten2[i] += M2(i,j) * ten[j];
          }
        }
      }
    }else if(step<=n_mem){ //step with expansion
      int last_back_block_size=1;
      for(int i=0; i<step-2; i++)last_back_block_size*=Ngrps2;

      for(int i=0; i<NL; i++){
        int gi=grp2[i];
        for(int j=0; j<NL; j++){ 
          int gj=grp2[j];
          for(int b=0; b<last_back_block_size; b++){
            ten2[(i*Ngrps2+gj)*last_back_block_size+b] +=
              IF.ten[step-1][(gi*Ngrps2+gj)*last_back_block_size+b] * 
              prop.M(i,j) * ten[j*last_back_block_size+b];
          }
        }
      }
      if(use_symmetric_Trotter){
        ten.swap(ten2);
        ten2.fill(0.);
        for(int i=0; i<NL; i++){
          for(int j=0; j<NL; j++){
            for(int b=0; b<last_back_block_size*Ngrps2; b++){
             ten2[i*last_back_block_size*Ngrps2+b] += 
                  M2(i,j) * ten[j*last_back_block_size*Ngrps2+b];
            }
          }
        }
      }
    }else{ //propagate and contract last index   
      int minus_two_block_size=1;
      for(int i=0; i<n_mem-2; i++)minus_two_block_size*=Ngrps2;

      for(int i=0; i<NL; i++){
        int gi=grp2[i];
        for(int j=0; j<NL; j++){ 
          int gj=grp2[j];
          for(int b=0; b<minus_two_block_size; b++){
            for(int l=0; l<Ngrps2; l++){ 
              ten2[(i*Ngrps2+gj)*minus_two_block_size+b] +=
                IF.ten[n_mem][((gi*Ngrps2+gj)*minus_two_block_size+b)*Ngrps2+l] * 
                prop.M(i,j) * ten[(j*minus_two_block_size+b)*Ngrps2+l];
            }
          }
        }
      }

      if(use_symmetric_Trotter){
        ten.swap(ten2);
        ten2.fill(0.);
        for(int i=0; i<NL; i++){
          for(int j=0; j<NL; j++){
            for(int b=0; b<minus_two_block_size*Ngrps2; b++){
             ten2[i*minus_two_block_size*Ngrps2+b] += 
                  M2(i,j) * ten[j*minus_two_block_size*Ngrps2+b];
            }
          }
        }
      }
    }

    ten.swap(ten2);

    update_rho(step);
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

    ten.fill(0.);
    for(int i=0; i<N; i++){
      for(int j=0; j<N; j++){
        ten[i*N+j]=rho(i,j);
      }
    }
  }

}
