#include "MPO_state.hpp"

namespace ACE{
namespace MPO_state{  //features with clear interpretation of MPO as state

std::vector<Eigen::VectorXcd> get_cap_d2(const MPO & mpo){
  if(mpo.size()<1)return std::vector<Eigen::VectorXcd>();
  std::vector<Eigen::VectorXcd> cap(mpo.size());
  cap[0]=Eigen::VectorXcd::Ones(mpo[0].dim_d2);
  for(int l=1; l<mpo.size(); l++){
    cap[l]=Eigen::VectorXcd::Zero(mpo[l-1].dim_d4);
    int N=sqrt(mpo[l-1].dim_d3);
    for(int d1=0; d1<mpo[l-1].dim_d1; d1++){
      for(int d2=0; d2<mpo[l-1].dim_d2; d2++){
        for(int i=0; i<N; i++){
          for(int d4=0; d4<mpo[l-1].dim_d4; d4++){
            cap[l](d4)+=cap[l-1](d2)*mpo[l-1](d1, d2, i*N+i, d4);
          }
        } 
      }
    }
  }
  return cap;
}
std::vector<Eigen::VectorXcd> get_cap_d4(const MPO & mpo){
  if(mpo.size()<1)return std::vector<Eigen::VectorXcd>();
  std::vector<Eigen::VectorXcd> cap(mpo.size());
  cap.back()=Eigen::VectorXcd::Ones(mpo.back().dim_d4);
  for(int l=(int)mpo.size()-2; l>=0; l--){
    cap[l]=Eigen::VectorXcd::Zero(mpo[l+1].dim_d2);
    int N=sqrt(mpo[l+1].dim_d3);
    for(int d1=0; d1<mpo[l+1].dim_d1; d1++){
      for(int d2=0; d2<mpo[l+1].dim_d2; d2++){
        for(int i=0; i<N; i++){
          for(int d4=0; d4<mpo[l+1].dim_d4; d4++){
            cap[l](d2)+=cap[l+1](d4)*mpo[l+1](d1, d2, i*N+i, d4);
          }
        } 
      }
    }
  }
  return cap;
}
std::complex<double> evaluate(const MPO & mpo, const std::vector<int> &list){
  if(mpo.size()!=list.size()){
    std::cerr<<"MPO_state::evaluate: mpo.size()!=list.size()!"<<std::endl;
    exit(1);
  }
  if(mpo.size()<1)return 0;
  if(mpo[0].dim_d2!=1){
    std::cerr<<"MPO_state::evaluate: mpo[0].dim_d2!=1!"<<std::endl;
    exit(1);
  }
  Eigen::VectorXcd cap=Eigen::VectorXcd::Ones(1);
  for(int l=0; l<mpo.size(); l++){
    if(mpo[l].dim_d1!=1){
      std::cerr<<"MPO_state::evaluate: mpo["<<l<<"].dim_d1!=1!"<<std::endl;
      exit(1);
    }
    if(mpo[l].dim_d2!=cap.rows()){
      std::cerr<<"MPO_state::evaluate: mpo["<<l<<"].dim_d2!=cap.rows()!"<<std::endl;
      exit(1);
    }
    if(list[l]>=mpo[l].dim_d3){
      std::cerr<<"MPO_state::evaluate: list["<<l<<"]>=mpo["<<l<<"].dim_d3!"<<std::endl;
      exit(1);
    }
    Eigen::VectorXcd cap2=Eigen::VectorXcd::Zero(mpo[l].dim_d4);
    for(int d2=0; d2<mpo[l].dim_d2; d2++){
      for(int d4=0; d4<mpo[l].dim_d4; d4++){
        cap2(d4)+=cap(d2)*mpo[l](0, d2, list[l], d4);
      }
    }
    cap2.swap(cap);
    if(mpo.back().dim_d4!=1){
      std::cerr<<"MPO_state::evaluate: mpo.back().dim_d4!=1!"<<std::endl;
      exit(1);
    }
  }
  return cap(0);
}

MPO from_SingleExcitationManifold(Eigen::MatrixXcd &rho, std::vector<int> dims){
//std::cout<<"TEST: rho="<<std::endl<<rho<<std::endl;
  int N=rho.rows();
  if(dims.size()==0){
    dims=std::vector<int>(N,2);
  }else if(dims.size()!=N){
    std::cerr<<"MPO_state::from_SingleExcitationManifold: dims.size()!=N!"<<std::endl;
    throw DummyException();
  }
  //initialize:
  MPO mpo(N);
  for(int r=0; r<N; r++){   
    if(dims[r]<2){
      std::cerr<<"MPO_state::from_SingleExcitationManifold: dims["<<r<<"<2!"<<std::endl;
      throw DummyException();
    }
    mpo[r]=FourLeg(1,1+N*N,dims[r]*dims[r],1+N*N); 
    mpo[r].set_zero();
  }  

  //elements:
  for(int i=0; i<N; i++){  
    for(int j=0; j<N; j++){  
      for(int l=0; l<N; l++){  
        if(l!=i && l!=j){ 
          mpo[l](0, 1+i*N+j, 0, 1+i*N+j)=1.;  //pipe
        }
        if(i==j && l==i){ //diag
          mpo[l](0, 1+i*N+j, 3 ,1+i*N+j)=rho(i,i); 
          continue; 
        } 
        if(l==i){ mpo[l](0,1+i*N+j,1*dims[l]+0,1+i*N+j)=sqrt(rho(i,j)); } //offdiag1
        if(l==j){ mpo[l](0,1+i*N+j,0*dims[l]+1,1+i*N+j)=sqrt(rho(i,j)); } //offdiag2
      }
    }
  } 
  Eigen::VectorXcd cap=Eigen::VectorXcd::Ones(1+N*N);
  mpo[0].multiply_d2(cap);
  mpo.back().multiply_d4(cap);

  return mpo;
}

std::complex<double> dot(const MPO & mpo1, const MPO & mpo2){
  std::complex<double> ret=0;
  if(mpo1.size()!=mpo2.size()){
    std::cerr<<"MPO_state::dot: mpo1.size()!=mpo2.size()!"<<std::endl;
    throw DummyException();
  }
  if(mpo1.size()<1)return 0;

  if(mpo1[0].dim_d2!=1){
    std::cerr<<"MPO_state::dot: mpo1[0].dim_d2!=1!"<<std::endl;
    throw DummyException();
  }   
  if(mpo2[0].dim_d2!=1){
    std::cerr<<"MPO_state::dot: mpo2[0].dim_d2!=1!"<<std::endl;
    throw DummyException();
  }  
  if(mpo1.back().dim_d4!=1){
    std::cerr<<"MPO_state::dot: mpo1.back().dim_d4!=1!"<<std::endl;
    throw DummyException();
  }   
  if(mpo2.back().dim_d4!=1){
    std::cerr<<"MPO_state::dot: mpo2.back().dim_d4!=1!"<<std::endl;
    throw DummyException();
  }
  Eigen::MatrixXcd cap(1,1); cap(0,0)=1.;
  for(int r=0; r<(int)mpo1.size(); r++){
    if(mpo1[r].dim_d1!=1){
      std::cerr<<"MPO_state::dot: mpo1[r].dim_d1!=1!"<<std::endl;
      throw DummyException();
    }   
    if(mpo2[r].dim_d1!=1){
      std::cerr<<"MPO_state::dot: mpo2[r].dim_d1!=1!"<<std::endl;
      throw DummyException();
    }    
    if(mpo1[r].dim_d3!=mpo2[r].dim_d3){
      std::cerr<<"MPO_state::dot: mpo1[r].dim_d3!=mpo2[r].dim_d3!"<<std::endl;
      throw DummyException();
    }      
    if(mpo1[r].dim_d2!=cap.rows()){
      std::cerr<<"MPO_state::dot: mpo1[r].dim_d2!=cap.rows()!"<<std::endl;
      throw DummyException();
    }    
    if(mpo2[r].dim_d2!=cap.cols()){
      std::cerr<<"MPO_state::dot: mpo2[r].dim_d2!=cap.cols()!"<<std::endl;
      throw DummyException();
    }    
    Eigen::MatrixXcd tmp=Eigen::MatrixXcd::Zero(mpo1[r].dim_d4*mpo1[r].dim_d3, cap.cols()); 
    //multiply cap with adjoint mpo1
    int N=sqrt(mpo1[r].dim_d3);
    for(int c=0; c<cap.cols(); c++){
      for(int d4=0; d4<mpo1[r].dim_d4; d4++){
        for(int i=0; i<N; i++){
          for(int j=0; j<N; j++){
            for(int p=0; p<cap.rows(); p++){
              tmp(d4*mpo1[r].dim_d3+(i*N+j),c)+=cap(p,c)*std::conj(mpo1[r](0,p,j*N+i,d4));
            }
          }
        }
      }
    }
    //multiply cap with mpo2
    cap=Eigen::MatrixXcd::Zero(mpo1[r].dim_d4, mpo2[r].dim_d4);
    for(int d42=0; d42<mpo2[r].dim_d4; d42++){
      for(int d41=0; d41<mpo1[r].dim_d4; d41++){
        for(int i=0; i<N; i++){
          for(int j=0; j<N; j++){
            for(int c=0; c<mpo2[r].dim_d2; c++){
              cap(d41,d42)+=tmp(d41*mpo1[r].dim_d3+(i*N+j),c)*mpo2[r](0,c,i*N+j,d42);
            }
          }
        }
      }
    }
  }
  return cap(0,0);
}

}//namespace MPO_state
}//namespace ACE
