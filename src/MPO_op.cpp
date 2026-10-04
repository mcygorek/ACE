#include "MPO_op.hpp"

namespace ACE{
namespace MPO_op{  //features with clear interpretation of MPO as operator

/*
  void add_coeff(MPO &mpo, const std::vector<int> & list_d1, const std::vector<int> & list_d3, const std::complex<double> & value){
    if(list_d1.size()!=list_d3.size()){
      std::cerr<<"MPO_op::add_coeff: list_d1.size()!=list_d3.size()!"<<std::endl;
      throw DummyException();
    }
    if(list_d1.size()!=mpo.size()){
      std::cerr<<"MPO_op::add_coeff: list_d1.size()!=mpo.size()!"<<std::endl;
      throw DummyException();
    } 
    if(mpo.size()<1){ return; }               
    if(mpo[0].dim_d2=!1){
      std::cerr<<"MPO_op::add_coeff: mpo[0].dim_d2=!1!"<<std::endl;
      throw DummyException();
    }                
    if(mpo.back().dim_d4=!1){
      std::cerr<<"MPO_op::add_coeff: mpo.back().dim_d4=!1!"<<std::endl;
      throw DummyException();
    }  
    if(mpo.size()==1){
      mpo[0](list_d1[0], 0, list_d3[0], 0)+=value;
      return;
    }
    std::complex<double> part_val=std::pow(value, 1./mpo.size());
    for(int l=0; l<mpo.size(); l++){
      if(list_d1[l]>=mpo[l].dim_d1){
        std::cerr<<"MPO_op::add_coeff: list_d1[l]>=mpo[l].dim_d1)!"<<std::endl;
        throw DummyException();
      }
      if(list_d3[l]>=mpo[l].dim_d3){
        std::cerr<<"MPO_op::add_coeff: list_d3[l]>=mpo[l].dim_d3)!"<<std::endl;
        throw DummyException();
      }
      { 
        FourLeg tmp; 
        tmp.resize_copy(mpo[l].dim_d1, mpo[l].dim_d2+1, mpo[l].dim_d3, mpo[l].dim_d4+1, mpo[l]);
        tmp.swap(mpo[l]);
      }

      mpo[l](list_d1[l], mpo[l].dim_d2-1, list_d3[l], mpo[l].dim_d4-1)=part_val;

      if(l==0){ mpo[l].multiply_d2(Eigen::VectorXcd::Ones(mpo[l].dim_d2)); }
      if(l==mpo.size()-1){ mpo[l].multiply_d4(Eigen::VectorXcd::Ones(mpo[l].dim_d4)); }
    }

  }
*/
std::complex<double> evaluate(const MPO & mpo, const std::vector<int> &list_d1, const std::vector<int> &list_d3){
  if(list_d3.size()!=list_d1.size()){
    std::cerr<<"MPO_op::evaluate: list_d3.size()!=list_d1.size()!"<<std::endl;
    throw DummyException();
  }
  if(mpo.size()!=list_d1.size()){
    std::cerr<<"MPO_op::evaluate: mpo.size()!=list_d1.size()!"<<std::endl;
    throw DummyException();
  }
  if(mpo.size()<1)return 0;
  if(mpo[0].dim_d2!=1){
    std::cerr<<"MPO_op::evaluate: mpo[0].dim_d2!=1!"<<std::endl;
    throw DummyException();
  }
  Eigen::VectorXcd cap=Eigen::VectorXcd::Ones(1);
  for(int l=0; l<mpo.size(); l++){
    if(mpo[l].dim_d2!=cap.rows()){
      std::cerr<<"MPO_op::evaluate: mpo["<<l<<"].dim_d2!=cap.rows()!"<<std::endl;
      throw DummyException();
    }
    if(list_d1[l]>=mpo[l].dim_d1){
      std::cerr<<"MPO_op::evaluate: list_d1["<<l<<"]>=mpo["<<l<<"].dim_d1!"<<std::endl;
      throw DummyException();
    }
    if(list_d3[l]>=mpo[l].dim_d3){
      std::cerr<<"MPO_op::evaluate: list_d3["<<l<<"]>=mpo["<<l<<"].dim_d3!"<<std::endl;
      throw DummyException();
    }
    Eigen::VectorXcd cap2=Eigen::VectorXcd::Zero(mpo[l].dim_d4);
    for(int d2=0; d2<mpo[l].dim_d2; d2++){
      for(int d4=0; d4<mpo[l].dim_d4; d4++){
        cap2(d4)+=cap(d2)*mpo[l](list_d1[l], d2, list_d3[l], d4);
      }
    }
    cap2.swap(cap);
    if(mpo.back().dim_d4!=1){
      std::cerr<<"MPO_op::evaluate: mpo.back().dim_d4!=1!"<<std::endl;
      throw DummyException();
    }
  }
  return cap(0);
}
std::complex<double> evaluate(const MPO & mpo, const std::vector<Eigen::VectorXcd> &closure_d1, const std::vector<int> &list_d3){
  if(list_d3.size()!=closure_d1.size()){
    std::cerr<<"MPO_op::evaluate: list_d3.size()!=closure_d1.size()!"<<std::endl;
    throw DummyException();
  }
  if(mpo.size()!=closure_d1.size()){
    std::cerr<<"MPO_op::evaluate: mpo.size()!=closure_d1.size()!"<<std::endl;
    throw DummyException();
  }
  if(mpo.size()<1)return 0;
  if(mpo[0].dim_d2!=1){
    std::cerr<<"MPO_op::evaluate: mpo[0].dim_d2!=1!"<<std::endl;
    throw DummyException();
  }
  Eigen::VectorXcd cap=Eigen::VectorXcd::Ones(1);
  for(int l=0; l<mpo.size(); l++){
    if(mpo[l].dim_d2!=cap.rows()){
      std::cerr<<"MPO_op::evaluate: mpo["<<l<<"].dim_d2!=cap.rows()!"<<std::endl;
      throw DummyException();
    }
    if(closure_d1[l].rows()!=mpo[l].dim_d1){
      std::cerr<<"MPO_op::evaluate: closure_d1["<<l<<"].rows()!=mpo["<<l<<"].dim_d1 ("<<closure_d1[l].rows()<<" vs. "<<mpo[l].dim_d1<<")!"<<std::endl;
      throw DummyException();
    }
    if(list_d3[l]>=mpo[l].dim_d3){
      std::cerr<<"MPO_op::evaluate: list_d3["<<l<<"]>=mpo["<<l<<"].dim_d3!"<<std::endl;
      throw DummyException();
    }
    Eigen::VectorXcd cap2=Eigen::VectorXcd::Zero(mpo[l].dim_d4);
    for(int d4=0; d4<mpo[l].dim_d4; d4++){
      for(int d2=0; d2<mpo[l].dim_d2; d2++){
        for(int d1=0; d1<mpo[l].dim_d1; d1++){
          cap2(d4)+=closure_d1[l](d1)*cap(d2)*mpo[l](d1, d2, list_d3[l], d4);
        }
      }
    }
    cap2.swap(cap);
    if(mpo.back().dim_d4!=1){
      std::cerr<<"MPO_op::evaluate: mpo.back().dim_d4!=1!"<<std::endl;
      throw DummyException();
    }
  }
  return cap(0);
}


std::vector<Eigen::VectorXcd> get_cap_d2(const MPO & mpo, const std::vector<Eigen::VectorXcd> &closure_d1){
  if(mpo.size()<1)return std::vector<Eigen::VectorXcd>();
  if(closure_d1.size()!=mpo.size()){
    std::cerr<<"MPO_op::get_cap_d2: closure_d1.size()!=mpo.size()!"<<std::endl;
    throw DummyException();
  }
  std::vector<Eigen::VectorXcd> cap(mpo.size());
  cap[0]=Eigen::VectorXcd::Ones(mpo[0].dim_d2);
  for(int l=1; l<mpo.size(); l++){
    if(closure_d1[l-1].rows()!=mpo[l-1].dim_d1){
      std::cerr<<"MPO_op::get_cap_d2: closure_d1["<<l-1<<"].rows()!=mpo["<<l-1<<"].dim_d1 ("<<closure_d1[l-1].rows()<<" vs. "<<mpo[l-1].dim_d1<<")!"<<std::endl;
      throw DummyException();
    }
    cap[l]=Eigen::VectorXcd::Zero(mpo[l-1].dim_d4);
    int N=sqrt(mpo[l-1].dim_d3);
    for(int d1=0; d1<mpo[l-1].dim_d1; d1++){
      for(int d2=0; d2<mpo[l-1].dim_d2; d2++){
        for(int i=0; i<N; i++){
          for(int d4=0; d4<mpo[l-1].dim_d4; d4++){
            cap[l](d4)+=closure_d1[l-1](d1)*cap[l-1](d2)*mpo[l-1](d1, d2, i*N+i, d4);
          }
        } 
      }
    }
  }
  return cap;
}
std::vector<Eigen::VectorXcd> get_cap_d4(const MPO & mpo, const std::vector<Eigen::VectorXcd> &closure_d1){
  if(mpo.size()<1)return std::vector<Eigen::VectorXcd>();
  if(closure_d1.size()!=mpo.size()){
    std::cerr<<"MPO_op::get_cap_d4: closure_d1.size()!=mpo.size()!"<<std::endl;
    throw DummyException();
  }
  std::vector<Eigen::VectorXcd> cap(mpo.size());
  cap.back()=Eigen::VectorXcd::Ones(mpo.back().dim_d4);
  for(int l=(int)mpo.size()-2; l>=0; l--){
    if(closure_d1[l+1].rows()!=mpo[l+1].dim_d1){
      std::cerr<<"MPO_op::get_cap_d4: closure_d1[l+1].rows()!=mpo[l+1].dim_d1!"<<std::endl;
      throw DummyException();
    }
    cap[l]=Eigen::VectorXcd::Zero(mpo[l+1].dim_d2);
    int N=sqrt(mpo[l+1].dim_d3);
    for(int d1=0; d1<mpo[l+1].dim_d1; d1++){
      for(int d2=0; d2<mpo[l+1].dim_d2; d2++){
        for(int i=0; i<N; i++){
          for(int d4=0; d4<mpo[l+1].dim_d4; d4++){
            cap[l](d2)+=closure_d1[l+1](d1)*cap[l+1](d4)*mpo[l+1](d1, d2, i*N+i, d4);
          }
        } 
      }
    }
  }
  return cap;
}


void join_d3_d1(MPO & mpo, const MPO & other){
  if(mpo.size()!=other.size()){
    std::cerr<<"MPO_op::join_d3_d1: mpo.size()!=other.size()!"<<std::endl;
    throw DummyException();
  }
  if(mpo.size()<1){return;}
  for(int l=0; l<(int)mpo.size(); l++){
    mpo[l].join_d3_d1(other[l]);
  }
}
void join_d3_d1_compress_d2_d4(MPO & mpo, const MPO & other, const TruncatedSVD &trunc, int verbosity){
  if(mpo.size()!=other.size()){
    std::cerr<<"MPO_op::join_d3_d1_compress_d2_d4: mpo.size()!=other.size()!"<<std::endl;
    throw DummyException();
  }
  if(mpo.size()<1){return;}

  Eigen::MatrixXcd P=Eigen::MatrixXcd::Identity(mpo[0].dim_d2*other[0].dim_d2,mpo[0].dim_d2*other[0].dim_d2);
  for(int l=0; l<(int)mpo.size()-1; l++){
    mpo[l].join_d3_d1_multiply_d2(other[l], P);
    //mpo[l].join_d3_d1(other[l]);
    //mpo[l].multiply_d2(P);
    Eigen::MatrixXcd A = mpo[l].get_Matrix_d1d2d3_d4();
    TruncatedSVD_Return ret = trunc.compress(A);
//    std::cout<<"["<<l<<"]: singular values: "<<ret.sigma.transpose()<<std::endl;
    double keep=ret.sigma(0);
    if(trunc.keep!=0.)keep=trunc.keep;
#ifdef PRINT_KEEP
std::cout<<"keep: "<<keep<<std::endl;
#endif
    mpo[l].set_from_Matrix_d1d2d3_d4( ret.U*keep, mpo[l].dim_d2, mpo[l].dim_d3);
    P=((ret.sigma/keep).asDiagonal()*ret.Vdagger).transpose();
    if(verbosity>0){std::cout<<"forward sweep of site "<<l<<": "<<A.rows()<<","<<A.cols()<<" -> "<<P.cols()<<std::endl;}
  }
  mpo.back().join_d3_d1_multiply_d2(other.back(), P);
  //mpo.back().join_d3_d1(other.back());
  //mpo.back().multiply_d2(P);
}

void add_coeff(MPO & mpo, const std::vector<int> &list_d1, const std::vector<int> &list_d3, const std::complex<double> &value){
  if(mpo.size()!=list_d1.size()){
    std::cerr<<"add_coeff: mpo.size()!=list_d1.size()!"<<std::endl;
    throw DummyException();
  }
  if(mpo.size()!=list_d3.size()){
    std::cerr<<"add_coeff: mpo.size()!=list_d3.size()!"<<std::endl;
    throw DummyException();
  }
  /*
  if(mpo[0].dim_d2!=1){
    std::cerr<<"add_coeff: mpo[0].dim_d2!=1!"<<std::endl;
    throw DummyException();
  }
  if(mpo.back().dim_d4!=1){
    std::cerr<<"add_coeff: mpo.back().dim_d4!=1!"<<std::endl;
    throw DummyException();
  }
  */
  if(mpo.size()==1){
    if(list_d1[0]<0||list_d1[0]>=mpo[0].dim_d1){
      std::cerr<<"add_coeff: list_d1[0]<0||list_d1[0]>=mpo[0].dim_d1!"<<std::endl;
      throw DummyException();
    }
    if(list_d3[0]<0||list_d3[0]>=mpo[0].dim_d3){
      std::cerr<<"add_coeff: list_d3[0]<0||list_d3[0]>=mpo[0].dim_d3!"<<std::endl;
      throw DummyException();
    }
    mpo[0](list_d1[0], 0, list_d3[0], 0)+=value;
    return;
  }

  // adding a "pipe" through the MPO (increasing dimension d2 and d4)
  std::complex<double> part_value=std::pow(value, 1./((int)mpo.size()));
  for(int r=0; r<(int)mpo.size(); r++){
    if(list_d1[r]<0||list_d1[r]>=mpo[r].dim_d1){ 
      std::cerr<<"add_coeff: list_d1["<<r<<"]<0||list_d1["<<r<<"]>=mpo["<<r<<"].dim_d1!"<<std::endl;
      throw DummyException();
    }
    if(list_d3[r]<0||list_d3[r]>=mpo[r].dim_d3){ 
      std::cerr<<"add_coeff: list_d3["<<r<<"]<0||list_d3["<<r<<"]>=mpo["<<r<<"].dim_d3!"<<std::endl;
      throw DummyException();
    }

    FourLeg tmp;
    if(r==0){
      tmp.resize_copy(mpo[0].dim_d1, mpo[0].dim_d2, mpo[0].dim_d3, mpo[0].dim_d4+1, mpo[0]);
    }else if(r==(int)mpo.size()-1){
      tmp.resize_copy(mpo[r].dim_d1, mpo[r].dim_d2+1, mpo[r].dim_d3, mpo[r].dim_d4, mpo[r]);
    }else{
      tmp.resize_copy(mpo[r].dim_d1, mpo[r].dim_d2+1, mpo[r].dim_d3, mpo[r].dim_d4+1, mpo[r]);
    }
    tmp(list_d1[r], tmp.dim_d2-1, list_d3[r], tmp.dim_d4-1)=part_value;
    mpo[r].swap(tmp);
  }
}
MPO H_to_L_forward(const MPO &mpo){
// expand mpo otimes Id_2
  MPO mpo_forward(mpo.size());
  for(int l=0; l<mpo.size(); l++){
    mpo_forward[l]=FourLeg::Zero(mpo[l].dim_d1*mpo[l].dim_d1, mpo[l].dim_d2, mpo[l].dim_d3*mpo[l].dim_d1, mpo[l].dim_d4);
    for(int d1=0; d1<mpo[l].dim_d1; d1++){
      for(int d2=0; d2<mpo[l].dim_d2; d2++){
        for(int d3=0; d3<mpo[l].dim_d3; d3++){
          for(int d4=0; d4<mpo[l].dim_d4; d4++){
            for(int i=0; i<mpo[l].dim_d1; i++){
              mpo_forward[l](d1*mpo[l].dim_d1+i, d2, d3*mpo[l].dim_d1+i, d4) =
               mpo[l](d1, d2, d3, d4);
            }
          }
        }
      }
    }
  }
  return mpo_forward;
}
MPO H_to_L_backward(const MPO &mpo){
// expand Id_2 otimes mpo
  MPO mpo_backward(mpo.size());
  for(int l=0; l<mpo.size(); l++){
    mpo_backward[l]=FourLeg::Zero(mpo[l].dim_d1*mpo[l].dim_d1, mpo[l].dim_d2, mpo[l].dim_d3*mpo[l].dim_d1, mpo[l].dim_d4);
    for(int d1=0; d1<mpo[l].dim_d1; d1++){
      for(int d2=0; d2<mpo[l].dim_d2; d2++){
        for(int d3=0; d3<mpo[l].dim_d3; d3++){
          for(int d4=0; d4<mpo[l].dim_d4; d4++){
            for(int i=0; i<mpo[l].dim_d1; i++){
              mpo_backward[l](i*mpo[l].dim_d1+d1, d2, i*mpo[l].dim_d3+d3, d4) =
               std::conj(mpo[l](d1, d2, d3, d4));
            }
          }
        }
      }
    }
  }
  return mpo_backward;
}
FourLeg join_H_to_L(const FourLeg & f1, const FourLeg & f2, bool conjugate_first){
  FourLeg f(f1.dim_d1*f2.dim_d1, f1.dim_d2*f2.dim_d2, f1.dim_d3*f2.dim_d3, f1.dim_d4*f2.dim_d4);
  if(conjugate_first){
    for(int d11=0; d11<f1.dim_d1; d11++){
     for(int d12=0; d12<f2.dim_d1; d12++){
      for(int d21=0; d21<f1.dim_d2; d21++){
       for(int d22=0; d22<f2.dim_d2; d22++){
        for(int d31=0; d31<f1.dim_d3; d31++){
         for(int d32=0; d32<f2.dim_d3; d32++){
          for(int d41=0; d41<f1.dim_d4; d41++){
           for(int d42=0; d42<f2.dim_d4; d42++){
             f(d11*f2.dim_d1+d12, d21*f2.dim_d2+d22,
                  d31*f2.dim_d3+d32, d41*f2.dim_d4+d42) =
                      std::conj(f1(d11,d21,d31,d41))*f2(d12,d22,d32,d42);
    }}}}}}}} 
  }else{ //conjugate second
    for(int d11=0; d11<f1.dim_d1; d11++){
     for(int d12=0; d12<f2.dim_d1; d12++){
      for(int d21=0; d21<f1.dim_d2; d21++){
       for(int d22=0; d22<f2.dim_d2; d22++){
        for(int d31=0; d31<f1.dim_d3; d31++){
         for(int d32=0; d32<f2.dim_d3; d32++){
          for(int d41=0; d41<f1.dim_d4; d41++){
           for(int d42=0; d42<f2.dim_d4; d42++){
             f(d11*f2.dim_d1+d12, d21*f2.dim_d2+d22,
                  d31*f2.dim_d3+d32, d41*f2.dim_d4+d42) =
                      f1(d11,d21,d31,d41)*std::conj(f2(d12,d22,d32,d42));
    }}}}}}}} 
  }
  return f;
}
MPO join_H_to_L(const MPO & mpo1, const MPO & mpo2, bool conjugate_first){
  // square propagator on Hilbert space to get propagator on Liouville space
  if(mpo1.size()!=mpo2.size()){
    std::cerr<<"MPO_op::join_H_to_L: mpo1.size()!=mpo2.size()!"<<std::endl;
    throw DummyException();
  }
  MPO mpo(mpo1.size());
  for(int l=0; l<mpo1.size(); l++){
    mpo[l]=join_H_to_L(mpo1[l], mpo2[l], conjugate_first);
  }
  return mpo;
}

MPO SingleExcitationManifold_to_MPO_H(const Eigen::MatrixXcd & H, const TruncatedSVD & trunc, std::vector<int> dims){
  int N=H.rows();
  if(dims.size()==0){
    dims=std::vector<int>(N,2);
  }else if(dims.size()!=N){
    std::cerr<<"SingleExcitationManifold_to_MPO_H: dims.size()!=N!"<<std::endl;
    throw DummyException();
  }
  //initialize:
  MPO mpo(N);
  for(int r=0; r<N; r++){   
    mpo[r]=FourLeg(dims[r],1+N*N,dims[r],1+N*N); 
    mpo[r].set_zero();
  }  
  //map H_ij to entry of d1=0,..,1_i,..0, d3=0,..,1_j,..,0
  // by adding a "pipe" through the MPO (increasing dimension d2 and d4)

  //TODO: per-site dimensions other than 2 NOT IMPLEMENTED YET!!!

  //elements:
  for(int i=0; i<N; i++){  
    for(int j=0; j<N; j++){  
      for(int l=0; l<N; l++){
        if(l!=i && l!=j){ mpo[l](0,1+i*N+j,0,1+i*N+j)=1.; } //pipe
        if(i==j && l==i){ mpo[l](1,1+i*N+j,1,1+i*N+j)=H(i,j); continue; } //diag
        if(l==j){ mpo[l](1,1+i*N+j,0,1+i*N+j)=sqrt(H(i,j)); } //offdiag1
        if(l==i){ mpo[l](0,1+i*N+j,1,1+i*N+j)=sqrt(H(i,j)); } //offdiag2
      }
    }
  } 
  Eigen::VectorXcd cap=Eigen::VectorXcd::Ones(1+N*N);
  mpo[0].multiply_d2(cap);
  mpo.back().multiply_d4(cap);

  //now compress
  mpo.sweep_d2_d4(trunc);
  mpo.sweep_d4_d2(trunc);
  return mpo;
}
std::complex<double> SingleExcitationManifold_evaluate_MPO(const MPO & mpo, int i, int j){
  Eigen::VectorXcd v=Eigen::VectorXcd::Ones(1);
  for(int l=0; l<mpo.size(); l++){
    int i1=0; 
    int i3=0;
    if(l==i){i1=1;}
    if(l==j){i3=1;}
    Eigen::VectorXcd v2=Eigen::VectorXcd::Zero(mpo[l].dim_d4);
    for(int d2=0; d2<mpo[l].dim_d2; d2++){
      for(int d4=0; d4<mpo[l].dim_d4; d4++){
         v2(d4) += mpo[l](i1, d2, i3, d4) * v(d2);     
      }
    }
    v.swap(v2);
  }
  return v(0);
}
double compare_SingeExcitationManifold_MPO(const MPO &mpo, const Eigen::MatrixXcd &M){
  if(mpo.size()!=M.rows()||mpo.size()!=M.cols()){
    std::cerr<<"MPO_op::compare_SingeExcitationManifold_MPO: mpo.size()!=M.rows()||mpo.size()!=M.cols()!"<<std::endl;
    throw DummyException();
  }
  Eigen::MatrixXcd Mrec(M.rows(), M.cols());
  for(int i=0; i<M.rows(); i++){
    for(int j=0; j<M.cols(); j++){
      Mrec(i,j)=MPO_op::SingleExcitationManifold_evaluate_MPO(mpo, i, j);
    }
  }
  return (Mrec-M).norm();
}
 
MPO Identity(const std::vector<int> & dims2){
  MPO mpo(dims2.size());
  for(size_t r=0; r<dims2.size(); r++){
    mpo[r]=FourLeg::Identity_d1_d3(dims2[r]);
  }
  return mpo;
}
MPO Zero(const std::vector<int> & dims2){
  MPO mpo(dims2.size());
  for(size_t r=0; r<dims2.size(); r++){
    mpo[r]=FourLeg::Zero(dims2[r],1,dims2[r],1);
  }
  return mpo;
}

MPO OneBody(const std::vector<int> & dims, int r1, const Eigen::MatrixXcd &op1){

  if(r1<0||r1>=dims.size()){
    std::cerr<<"MPO_op::OneBody: r1<0||r1>=dims.size()!"<<std::endl;
    throw DummyException();
  }
  if(op1.rows()!=dims[r1] ){
    std::cerr<<"MPO_op::OneBody: op1.rows()!=dims[r1]!"<<std::endl;
    throw DummyException();
  }
  if(op1.cols()!=dims[r1] ){
    std::cerr<<"MPO_op::OneBody: op1.cols()!=dims[r1]!"<<std::endl;
    throw DummyException();
  }
  
  MPO mpo=MPO_op::Identity(dims);
  mpo[r1]=FourLeg::Operator_d1_d3(op1);
  return mpo;
}
MPO OneBody_forward(const std::vector<int> & dims2, int r1, const Eigen::MatrixXcd &op1){

  if(r1<0||r1>=dims2.size()){
    std::cerr<<"MPO_op::OneBody_forward: r1<0||r1>=dims2.size()!"<<std::endl;
    throw DummyException();
  }
  if(op1.rows()*op1.rows()!=dims2[r1] ){
    std::cerr<<"MPO_op::OneBody_forward: op1.rows()^2!=dims2[r1]!"<<std::endl;
    throw DummyException();
  }
  if(op1.cols()*op1.cols()!=dims2[r1] ){
    std::cerr<<"MPO_op::OneBody_forward: op1.cols()^2!=dims2[r1]!"<<std::endl;
    throw DummyException();
  }
  
  MPO mpo=MPO_op::Identity(dims2);
  mpo[r1]=FourLeg::Operator_forward_d1_d3(op1);
  return mpo;
}
MPO OneBody_backward(const std::vector<int> & dims2, int r1, const Eigen::MatrixXcd &op1){

  if(r1<0||r1>=dims2.size()){
    std::cerr<<"MPO_op::OneBody_backward: r1<0||r1>=dims2.size()!"<<std::endl;
    throw DummyException();
  } 
  if(op1.rows()*op1.rows()!=dims2[r1]){
    std::cerr<<"MPO_op::OneBody_backward: op1.rows()^2!=dims2[r1]!"<<std::endl;
    throw DummyException();
  }
  if(op1.cols()*op1.cols()!=dims2[r1]){
    std::cerr<<"MPO_op::OneBody_backward: op1.cols()^2!=dims2[r1]!"<<std::endl;
    throw DummyException();
  }
  
  MPO mpo=MPO_op::Identity(dims2);
  mpo[r1]=FourLeg::Operator_backward_d1_d3(op1);
  return mpo;
}

MPO OneBody_forward_backward(const std::vector<int> & dims2, int r1, const Eigen::MatrixXcd &op1, const Eigen::MatrixXcd &op2){

  if(r1<0||r1>=dims2.size()){
    std::cerr<<"MPO_op::OneBody_forward_backward: r1<0||r1>=dims2.size()!"<<std::endl;
    throw DummyException();
  }
  if(op1.rows()*op1.rows()!=dims2[r1] ){
    std::cerr<<"MPO_op::OneBody_forward_backward: op1.rows()^2!=dims2[r1]!"<<std::endl;
    throw DummyException();
  }
  if(op1.cols()*op1.cols()!=dims2[r1] ){
    std::cerr<<"MPO_op::OneBody_forward_backward: op1.cols()^2!=dims2[r1]!"<<std::endl;
    throw DummyException();
  }
  if(op2.rows()!=op1.cols() || op2.cols()!=op1.rows()){
    std::cerr<<"MPO_op::OneBody_forward_backward: op2.rows()!=op1.cols() || op2.cols()!=op1.rows()!"<<std::endl;
    throw DummyException();
  }
  
  MPO mpo=MPO_op::Identity(dims2);
  mpo[r1]=FourLeg::Operator_forward_backward_d1_d3(op1, op2);
  return mpo;
}

MPO OneBody_Sum(const std::vector<Eigen::MatrixXcd> & list){
  std::vector<int> dims(list.size());
  for(size_t i=0; i<list.size(); i++){
    dims[i]=list[i].rows();
  }
  MPO mpo = Zero(dims);
  for(size_t i=0; i<list.size(); i++){
    MPO tmp = OneBody(dims, i, list[i]);
    add(mpo, tmp);
  }
  return mpo;
}


MPO TwoBody(const std::vector<int> & dims, int r1, int r2, 
                   const Eigen::MatrixXcd &op1, const Eigen::MatrixXcd &op2){

  if(r1<0||r1>=dims.size() || r2<0||r2>=dims.size()){
    std::cerr<<"MPO_op::TwoBody: r1<0||r1>=dims.size() || r2<0||r2>=dims.size()!"<<std::endl;
    throw DummyException();
  }
  if(r1==r2){
    std::cerr<<"MPO_op::TwoBody: r1==r2!"<<std::endl;
    throw DummyException();
  } 
  if(op1.rows()!=dims[r1] || op2.rows()!=dims[r2]){
    std::cerr<<"MPO_op::TwoBody: op1.rows()!=dims[r1] || op2.rows()!=dims[r2]!"<<std::endl;
    throw DummyException();
  }
  if(op1.cols()!=dims[r1] || op2.cols()!=dims[r2]){
    std::cerr<<"MPO_op::TwoBody: op1.cols()!=dims[r1] || op2.cols()!=dims[r2]!"<<std::endl;
    throw DummyException();
  }
  
  MPO mpo=MPO_op::Identity(dims);
  mpo[r1]=FourLeg::Operator_d1_d3(op1);
  mpo[r2]=FourLeg::Operator_d1_d3(op2);
  return mpo;
}
MPO TwoBody_forward(const std::vector<int> & dims2, int r1, int r2, 
                   const Eigen::MatrixXcd &op1, const Eigen::MatrixXcd &op2){

  if(r1<0||r1>=dims2.size() || r2<0||r2>=dims2.size()){
    std::cerr<<"MPO_op::TwoBody_forward: r1<0||r1>=dims2.size() || r2<0||r2>=dims2.size()!"<<std::endl;
    throw DummyException();
  }
  if(r1==r2){
    std::cerr<<"MPO_op::TwoBody_forward: r1==r2!"<<std::endl;
    throw DummyException();
  } 
  if(op1.rows()*op1.rows()!=dims2[r1] || op2.rows()*op2.rows()!=dims2[r2]){
    std::cerr<<"MPO_op::TwoBody_forward: op1.rows()^2!=dims2[r1] || op2.rows()^2!=dims2[r2]!"<<std::endl;
    throw DummyException();
  }
  if(op1.cols()*op1.cols()!=dims2[r1] || op2.cols()*op2.cols()!=dims2[r2]){
    std::cerr<<"MPO_op::TwoBody_forward: op1.cols()^2!=dims2[r1] || op2.cols()^2!=dims2[r2]!"<<std::endl;
    throw DummyException();
  }
  
  MPO mpo=MPO_op::Identity(dims2);
  mpo[r1]=FourLeg::Operator_forward_d1_d3(op1);
  mpo[r2]=FourLeg::Operator_forward_d1_d3(op2);
  return mpo;
}
MPO TwoBody_backward(const std::vector<int> & dims2, int r1, int r2, 
                   const Eigen::MatrixXcd &op1, const Eigen::MatrixXcd &op2){

  if(r1<0||r1>=dims2.size() || r2<0||r2>=dims2.size()){
    std::cerr<<"MPO_op::TwoBody_backward: r1<0||r1>=dims2.size() || r2<0||r2>=dims2.size()!"<<std::endl;
    throw DummyException();
  }
  if(r1==r2){
    std::cerr<<"MPO_op::TwoBody_backward: r1==r2!"<<std::endl;
    throw DummyException();
  } 
  if(op1.rows()*op1.rows()!=dims2[r1] || op2.rows()*op2.rows()!=dims2[r2]){
    std::cerr<<"MPO_op::TwoBody_backward: op1.rows()^2!=dims2[r1] || op2.rows()^2!=dims2[r2]!"<<std::endl;
    throw DummyException();
  }
  if(op1.cols()*op1.cols()!=dims2[r1] || op2.cols()*op2.cols()!=dims2[r2]){
    std::cerr<<"MPO_op::TwoBody_backward: op1.cols()^2!=dims2[r1] || op2.cols()^2!=dims2[r2]!"<<std::endl;
    throw DummyException();
  }
  
  MPO mpo=MPO_op::Identity(dims2);
  mpo[r1]=FourLeg::Operator_backward_d1_d3(op1);
  mpo[r2]=FourLeg::Operator_backward_d1_d3(op2);
  return mpo;
}
MPO TwoBody_forward_backward(const std::vector<int> & dims2, int r1, int r2, 
                   const Eigen::MatrixXcd &op1, const Eigen::MatrixXcd &op2){

  if(r1<0||r1>=dims2.size() || r2<0||r2>=dims2.size()){
    std::cerr<<"MPO_op::TwoBody_forward_backward: r1<0||r1>=dims2.size() || r2<0||r2>=dims2.size()!"<<std::endl;
    throw DummyException();
  }
  if(r1==r2){
    std::cerr<<"MPO_op::TwoBody_forward_backward: r1==r2!"<<std::endl;
    throw DummyException();
  } 
  if(op1.rows()*op1.rows()!=dims2[r1] || op2.rows()*op2.rows()!=dims2[r2]){
    std::cerr<<"MPO_op::TwoBody_forward_backward: op1.rows()^2!=dims2[r1] || op2.rows()^2!=dims2[r2]!"<<std::endl;
    throw DummyException();
  }
  if(op1.cols()*op1.cols()!=dims2[r1] || op2.cols()*op2.cols()!=dims2[r2]){
    std::cerr<<"MPO_op::TwoBody_forward_backward: op1.cols()^2!=dims2[r1] || op2.cols()^2!=dims2[r2]!"<<std::endl;
    throw DummyException();
  }
  
  MPO mpo=MPO_op::Identity(dims2);
  mpo[r1]=FourLeg::Operator_forward_d1_d3(op1);
  mpo[r2]=FourLeg::Operator_backward_d1_d3(op2);
  return mpo;
}

MPO Product(const std::vector<Eigen::MatrixXcd> & list){
  MPO mpo;
  mpo.resize(list.size());
  for(size_t i=0; i<list.size(); i++){
    mpo[i]=FourLeg::Operator_d1_d3(list[i]);
  }
  return mpo;
}


void add(MPO & mpo, const MPO & other){
  if(mpo.size()!=other.size()){
    std::cerr<<"MPO_op::add: mpo.size()!=other.size()!"<<std::endl;
    throw DummyException();
  }
  if(mpo.size()<1)return;
  if(mpo[0].dim_d2!=1){
    std::cerr<<"MPO_op::add: mpo[0].dim_d2="<<mpo[0].dim_d2<<"!=1!"<<std::endl;
    throw DummyException();
  }
  if(other[0].dim_d2!=1){
    std::cerr<<"MPO_op::add: other[0].dim_d2="<<other[0].dim_d2<<"!=1!"<<std::endl;
    throw DummyException();
  }
  if(mpo.back().dim_d4!=1 || other.back().dim_d4!=1){
    std::cerr<<"MPO_op::add: mpo.back().dim_d4!=1 || other.back().dim_d4!=1!"<<std::endl;
    throw DummyException();
  }
  for(size_t r=0; r<other.size(); r++){
    if(mpo[r].dim_d1!=other[r].dim_d1){
      std::cerr<<"MPO_op::add: mpo[r].dim_d1!=other[r].dim_d1!"<<std::endl;
      throw DummyException();
    }
    if(mpo[r].dim_d3!=other[r].dim_d3){
      std::cerr<<"MPO_op::add: mpo[r].dim_d3!=other[r].dim_d3!"<<std::endl;
      throw DummyException();
    }
    FourLeg f=FourLeg::Zero(mpo[r].dim_d1, mpo[r].dim_d2+other[r].dim_d2, mpo[r].dim_d3, mpo[r].dim_d4+other[r].dim_d4);
    for(int d1=0; d1<other[r].dim_d1; d1++){
      for(int d3=0; d3<other[r].dim_d3; d3++){
        for(int d2=0; d2<mpo[r].dim_d2; d2++){
          for(int d4=0; d4<mpo[r].dim_d4; d4++){
            f(d1, d2+other[r].dim_d2, d3, d4+other[r].dim_d4) = mpo[r](d1, d2, d3, d4);
          } 
        }
        for(int d2=0; d2<other[r].dim_d2; d2++){
          for(int d4=0; d4<other[r].dim_d4; d4++){
            f(d1, d2, d3, d4) = other[r](d1, d2, d3, d4);
          } 
        }
      } 
    }
    mpo[r].swap(f);
  } 
  mpo[0].multiply_d2(Eigen::VectorXcd::Ones(mpo[0].dim_d2));
  mpo.back().multiply_d4(Eigen::VectorXcd::Ones(mpo.back().dim_d4));
}

MPO exp_Taylor(const MPO & other, const std::complex<double> & scale, const TruncatedSVD & trunc, int N_Taylor){
  std::vector<int> dims2(other.size());
  for(size_t r=0; r<other.size(); r++){ 
    if(other[r].dim_d1!=other[r].dim_d3){
      std::cerr<<"MPO_op::exp: other[r].dim_d1!=other[r].dim_d3!"<<std::endl;
      throw DummyException();
    }
    dims2[r]=other[r].dim_d1;
  }
  MPO mpo = MPO_op::Identity(dims2);
  MPO tmp = other; 
  tmp.scale_first(scale);
  MPO_op::add(mpo, tmp);
  for(int n=2; n<=N_Taylor; n++){
    tmp.scale_first(scale/((double)n));
    MPO_op::join_d3_d1_compress_d2_d4(tmp, other, trunc);
    tmp.sweep_d4_d2(trunc);
    MPO_op::add(mpo, tmp);
    mpo.sweep_d2_d4(trunc);
    mpo.sweep_d4_d2(trunc);
  }
  return mpo;
}

}//namespace MPO_op
}//namespace ACE

