#include "ADM_MPS.hpp"
#include "Propagator.hpp"
#include "InfluenceFunctional_Vector.hpp"
#include "Tensor.hpp"
#include "LiouvilleTools.hpp"

namespace ACE{

  void AugmentedDensityMatrix_MPS::check_dimensions()const{
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

  void AugmentedDensityMatrix_MPS::print_status(std::ostream &os)const{
    os<<"ADM MPS max_dim: "<<ten.get_max_dim()<<std::endl;
  }

  void AugmentedDensityMatrix_MPS::update_rho(int step){
    int n_mem_eff=ten.get_rank();

    int NL=get_NL();
    int N=sqrt(NL);

//    int Ngrps2=get_Ngrps2();

    if(n_mem_eff<1){
      std::cerr<<"n_mem_eff must be larger than zero!"<<std::endl;
      exit(1);
    }
    if(n_mem_eff<2){
      rho=Eigen::MatrixXcd::Zero(N,N);
      for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
          rho(i,j)+=ten.a[0](i*N+j,0,0);
        }
      }
      return;
    }


    std::vector<std::complex<double> > A,B;
    A.resize(ten.a[n_mem_eff-1].dim_d1,0.);
    for(int i=0; i<ten.a[n_mem_eff-1].dim_i; i++){
      for(int d=0; d<ten.a[n_mem_eff-1].dim_d1; d++){
        A[d]+=ten.a[n_mem_eff-1](i, d, 0);
      }
    }
    
    for(int n=n_mem_eff-2; n>0; n--){
      B.clear();
      B.resize(ten.a[n].dim_d1, 0.);
      for(int d1=0; d1<ten.a[n].dim_d1; d1++){
        for(int i=0; i<ten.a[n].dim_i; i++){
          for(int d2=0; d2<ten.a[n].dim_d2; d2++){
            B[d1]+=ten.a[n](i,d1,d2)*A[d2];
          } 
        }
      }
      A.swap(B);
    }



    rho=Eigen::MatrixXcd::Zero(N,N);
    for(int i=0; i<N; i++){
      for(int j=0; j<N; j++){
        for(int d2=0; d2<ten.a[0].dim_d2; d2++){
          rho(i,j)+=ten.a[0](i*N+j,0,d2)*A[d2];
        }
      }
    }
  }

  void AugmentedDensityMatrix_MPS::propagate(Propagator &prop, 
           const InfluenceFunctional_Vector &IF, double t, double dt, 
           int step, RankCompressor *compressor, bool use_symmetric_Trotter){

    check_dimensions();
    int n_max=get_n_max();
    int n_mem_IF=IF.get_n_max();
std::cout<<"step="<<step<<" n_max="<<n_max<<" n_mem_IF="<<n_mem_IF<<std::endl;
    if(n_mem_IF<n_max+1){
      std::cerr<<"ADM_MPS::propagate: n_mem_IF<n_max+1!"<<std::endl;
      exit(1);
    }
 
    int NL=get_NL();
    int N=sqrt(NL);
    if(N!=IF.groups.sys_dim()){
      std::cerr<<"ADM::propagate: N!=IF.groups.sys_dim()!"<<std::endl;
      exit(1);
    }

    int Ngrps2=IF.lgroups.get_Ngrps();
    const std::vector<int> & grp2=IF.lgroups.grp;

    if(use_symmetric_Trotter){
      if(step==1){
        prop.update(t,dt/2.);
      }else{
        prop.update(t-dt/2.,dt);
      }
    }else{
      prop.update(t,dt);
    }

    MPS ten2;

    if(step<1){
      std::cerr<<"AugmentedDensityMatrix::propagate: step <1!"<<std::endl;
      exit(1);
    }else if(step == 1 || n_mem_IF<2){//first step: contract last index
      ten2.a.resize(1);
      ten2.a[0].resize(NL,1,1);
      ten2.a[0].fill(0.);
      for(int i=0; i<NL; i++){
        int gi=grp2[i];
        for(int j=0; j<NL; j++){
          ten2.a[0](i,0,0) += IF.b[0](gi,gi)* prop.M(i,j) * ten.a[0](j,0,0);
        }
      }
    }else if(step==n_mem_IF && n_mem_IF==2){ //M and truncation over same index
      ten2.a.resize(1);
      ten2.a[0].resize(NL,1,1);
      ten2.a[0].fill(0.);
      for(int i=0; i<NL; i++){
        int gi=grp2[i];
        for(int j=0; j<NL; j++){
          int gj=grp2[j];
          ten2.a[0](i,0,0) += IF.b[0](gi,gi) * IF.b[1](gi,gj) 
                              * prop.M(i,j) * ten.a[0](j,0,0);
        }
      }
    }else{

      ten2.a.resize(ten.a.size()+1);
      ten2.a[0].resize(NL,1,NL);
      ten2.a[0].fill(0.);
      for(int i=0; i<NL; i++){
        int gi=grp2[i];
/*
        for(int j=0; j<NL; j++){
          int gj=grp2[j];
          ten2.a[0](i, 0, j)=prop.M(i,j)*IF.b[0](gi,gi) ;
        }
*/
          ten2.a[0](i, 0, i)=IF.b[0](gi,gi) ;
      } 

      ten2.a[1].resize(Ngrps2, NL, Ngrps2*ten.a[0].dim_d2);
      ten2.a[1].fill(0.);
      for(int i=0; i<NL; i++){
        int gi=grp2[i];
        for(int j=0; j<NL; j++){ 
          int gj=grp2[j];
          for(int d=0; d<ten.a[0].dim_d2; d++){
            ten2.a[1](gi, j, d*Ngrps2+gj) +=
              IF.b[1](gj, gi) * prop.M(j,i)* ten.a[0](i, 0, d);
          }
        }
      }
      for(int n=2; n<(int)ten2.a.size(); n++){
        ten2.a[n].resize(Ngrps2, Ngrps2*ten.a[n-1].dim_d1, Ngrps2*ten.a[n-1].dim_d2);
        ten2.a[n].fill(0.);
        for(int k=0; k<Ngrps2; k++){ 
          for(int d1=0; d1<ten.a[n-1].dim_d1; d1++){
            for(int d2=0; d2<ten.a[n-1].dim_d2; d2++){
              for(int l=0; l<Ngrps2; l++){
                ten2.a[n](k, d1*Ngrps2+l, d2*Ngrps2+l ) +=
                  IF.b[n](l, k) * ten.a[n-1](k, d1, d2);
              }
            }
          }
        }
      }
      //last index: sum over delta:
      {
        MPS_Matrix &ref=ten2.a[ten2.a.size()-1];
        MPS_Matrix mmat(ref.dim_i, ref.dim_d1, 1);
        mmat.fill(0.);
        for(int i=0; i<ref.dim_i; i++){
          for(int d1=0; d1<ref.dim_d1; d1++){
            for(int d2=0; d2<ref.dim_d2; d2++){
              mmat(i,d1,0)+=ref(i,d1,d2);
            } 
          }
        }
        ref.swap(mmat);
      }

      if(step>n_max){ //contract last index
        MPS_Matrix &ref=ten2.a[ten2.a.size()-2];
        MPS_Matrix mmat(ref.dim_i, ref.dim_d1, 1);
        mmat.fill(0.);
        for(int i=0; i<ref.dim_i; i++){
          for(int d1=0; d1<ref.dim_d1; d1++){
            for(int d2=0; d2<ref.dim_d2; d2++){
              for(int i2=0; i2<ten2.a.back().dim_i; i2++){
                mmat(i,d1,0)+=ref(i,d1,d2)*ten2.a.back()(i2,d2,0);
              }
            }
          }
        }
        ten2.a[ten2.a.size()-2].swap(mmat);
      }
 
    }
 

    //swap ten with ten2
    if(step<=n_max){ 
      ten.swap(ten2);
    }else{
      for(size_t n=0; n<ten.a.size(); n++){
        ten.a[n].swap(ten2.a[n]);
      }
    }

    if(compressor!=NULL){
      if(compressor->debug){
        ten.check_consistency();
        std::cout<<"Inner indices before sweeps:"<<std::endl;
        for(int i=0; i<ten.get_rank(); i++){
          std::cout<<ten.a[i].dim_d1<<" ";
        }
        std::cout<<std::endl;
      }
 
      for(int i=0; i<(int)ten.a.size()-1; i++){
        compressor->sweep_block_low_to_high(i, ten);
      }
      for(int i=ten.a.size()-1; i>1; i--){
        compressor->sweep_block_high_to_low(i, ten);
      }

      if(compressor->debug){
        std::cout<<"Inner indices after sweeps:"<<std::endl;
        for(int i=0; i<ten.get_rank(); i++){
          std::cout<<ten.a[i].dim_d1<<" ";
        }
        std::cout<<std::endl;
      }
    }

    update_rho(step);
    if(use_symmetric_Trotter){
      prop.update(t+dt/2.,dt/2.);
      rho=L_Vector_to_H_Matrix( prop.M * H_Matrix_to_L_Vector(rho) );
    }
  }
  
  
  AugmentedDensityMatrix_MPS::AugmentedDensityMatrix_MPS(int n_max_, 
         int Ngrps, const Eigen::MatrixXcd &rho_)  : rho(rho_), n_max(n_max_){

    int N=rho_.rows();
    if(rho.rows()!=rho.cols() ){
      std::cerr<<"ADM: constructor: rho has wrong dimensions!"<<std::endl;
      exit(1);
    }

    if(n_max<0){
      std::cerr<<"n_max<0!"<<std::endl;
      exit(1);
    }
    ten.a.resize(1);
    ten.a[0].resize(N*N,1,1);
    for(int i=0; i<N; i++){
      for(int j=0; j<N; j++){
        ten.a[0](i*N+j,0,0)=rho(i,j);
      }
    }
  }
}
