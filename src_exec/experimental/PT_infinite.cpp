#include "ACE.hpp"
#include "ProcessTensorBuffer.hpp"
#include "ProcessTensorRepeat.hpp"

using namespace ACE;


Eigen::MatrixXcd PT_iTEBD_calc_R(const Eigen::VectorXcd & LambdaB, const Eigen::MatrixXcd &GLG2){
  int dim_d=LambdaB.size();
  int dim_i=GLG2.rows()/dim_d;
  Eigen::MatrixXcd R=Eigen::MatrixXcd::Zero(dim_d*dim_d, dim_d*dim_d);
  for(int i=0; i<dim_i; i++){
    for(int j=0; j<dim_i; j++){
      for(int d1=0; d1<dim_d; d1++){
        for(int d2=0; d2<dim_d; d2++){
          for(int d3=0; d3<dim_d; d3++){
            for(int d4=0; d4<dim_d; d4++){
              R(d1*dim_d+d2, d3*dim_d+d4)+=
GLG2(i*dim_d+d1, j*dim_d+d3)*LambdaB(d3)*std::conj(GLG2(i*dim_d+d2,j*dim_d+d4)*LambdaB(d4));
            }
          }
        }
      }
    }
  }
  return R;
}
Eigen::MatrixXcd PT_iTEBD_calc_L(const Eigen::VectorXcd & LambdaB, const Eigen::MatrixXcd &GLG2){
  int dim_d=LambdaB.size();
  int dim_i=GLG2.rows()/dim_d;
  Eigen::MatrixXcd L=Eigen::MatrixXcd::Zero(dim_d*dim_d, dim_d*dim_d);
  for(int i=0; i<dim_i; i++){
    for(int j=0; j<dim_i; j++){
      for(int d1=0; d1<dim_d; d1++){
        for(int d2=0; d2<dim_d; d2++){
          for(int d3=0; d3<dim_d; d3++){
            for(int d4=0; d4<dim_d; d4++){
              L(d1*dim_d+d2, d3*dim_d+d4)+=
LambdaB(d1)*GLG2(i*dim_d+d1, j*dim_d+d3)*std::conj(LambdaB(d2)*GLG2(i*dim_d+d2,j*dim_d+d4));
            }
          }
        }
      }
    }
  }
  return L;
}
/*
Eigen::MatrixXcd PT_iTEBD_right(const Eigen::VectorXcd & LambdaB, const Eigen::MatrixXcd &GLG2, int IMAX=100, double thr=1e-13){
    //For now: power iteration
  int dim_d=LambdaB.size();
  int dim_i=GLG2.rows()/dim_d;
std::cout<<"dim_i="<<dim_i<<" dim_d="<<dim_d<<std::endl;
  Eigen::MatrixXcd VR=Eigen::MatrixXcd::Ones(dim_d, dim_d);
  for(int iter=0; iter<IMAX; iter++){
std::cout<<"iter: "<<iter<<" VR.norm()="<<VR.norm()<<std::endl;
    Eigen::MatrixXcd tmp=Eigen::MatrixXcd::Zero(dim_i*dim_d, dim_i*dim_d);
    for(int i=0; i<dim_i; i++){
      for(int j=0; j<dim_i; j++){
        for(int d1=0; d1<dim_d; d1++){
          for(int d2=0; d2<dim_d; d2++){
            for(int d3=0; d3<dim_d; d3++){
              tmp(i*dim_d+d1, j*dim_d+d2)+=
                GLG2(i*dim_d+d1, j*dim_d+d3)*LambdaB(d3)*VR(d3, d2);
            }
          }
        }
      }
    }
std::cout<<"iter: "<<iter<<" tmp.norm()="<<tmp.norm()<<std::endl;
    Eigen::MatrixXcd VR2=Eigen::MatrixXcd::Zero(dim_d, dim_d);
    for(int i=0; i<dim_i; i++){
      for(int j=0; j<dim_i; j++){
        for(int d1=0; d1<dim_d; d1++){
          for(int d2=0; d2<dim_d; d2++){
            for(int d3=0; d3<dim_d; d3++){
              VR2(d1, d2)+=
tmp(i*dim_d+d1, j*dim_d+d3)*std::conj(GLG2(i*dim_d+d2, j*dim_d+d3)*LambdaB(d3));
            }
          }
        }
      }
    }
    double norm=VR2.norm();
    std::cout<<"norm="<<norm<<std::endl;
    VR2/=norm; 
    double diff=(VR-VR2).norm();
    VR=VR2;
    std::cout<<"diff("<<iter<<"): "<<diff<<std::endl;
    if(diff<thr){break};
  }
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solve(VR);
  Eigen::VectorXcd sqrt_D(dim_d);
  for(int d=0; d<dim_d; d++){  
    sqrt_D(d)=sqrt((std::complex<double>)solve.eigenvalues()(d));
  }
  Eigen::MatrixXcd X=solve.eigenvectors()*sqrt_D.asDiagonal();
  std::cout<<"(VR-X*X.adjoint()).norm()="<<(VR-X*X.adjoint()).norm()<<std::endl;
  std::cout<<"X="<<X<<std::endl;
  return X;
}
*/
std::pair<Eigen::MatrixXcd, Eigen::VectorXcd> PT_iTEBD_X(const Eigen::VectorXcd & LambdaB, const Eigen::MatrixXcd &GLG2){
  int dim_d=LambdaB.size();
  Eigen::MatrixXcd R=PT_iTEBD_calc_R(LambdaB, GLG2);
  //Largest eigenvalue -> matrix VR:
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solveR(R);
  int imax=0;
  for(int i=1; i<dim_d*dim_d; i++){
    if(std::abs(solveR.eigenvalues()(i))>std::abs(solveR.eigenvalues()(imax))){
      imax=i;
    }
  }
  Eigen::MatrixXcd VR(dim_d, dim_d);
  for(int d1=0; d1<dim_d; d1++){ 
    for(int d2=0; d2<dim_d; d2++){
      VR(d1,d2)=solveR.eigenvectors()(d1*dim_d+d2,imax);
    }
  } 
  //Decompose VR:
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solve(VR);
  Eigen::VectorXcd sqrt_D(dim_d);
  for(int d=0; d<dim_d; d++){  
    sqrt_D(d)=sqrt((std::complex<double>)solve.eigenvalues()(d));
  }
  return std::make_pair(solve.eigenvectors(), sqrt_D);
}

std::pair<Eigen::MatrixXcd, Eigen::VectorXcd> PT_iTEBD_Y(const Eigen::VectorXcd & LambdaB, const Eigen::MatrixXcd &GLG2){
  int dim_d=LambdaB.size();
  Eigen::MatrixXcd L=PT_iTEBD_calc_L(LambdaB, GLG2);
  //Largest eigenvalue -> matrix VR:
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solveL(L.transpose());
  int imax=0;
  for(int i=1; i<dim_d*dim_d; i++){
    if(std::abs(solveL.eigenvalues()(i))>std::abs(solveL.eigenvalues()(imax))){
      imax=i;
    }
  }
  Eigen::MatrixXcd VL(dim_d, dim_d);
  for(int d1=0; d1<dim_d; d1++){
    for(int d2=0; d2<dim_d; d2++){
      VL(d1,d2)=solveL.eigenvectors()(d1*dim_d+d2,imax);
    }
  }
  //Decompose VL:
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solve(VL);

  Eigen::VectorXcd sqrt_D(dim_d);
  for(int d=0; d<dim_d; d++){  
    sqrt_D(d)=sqrt((std::complex<double>)solve.eigenvalues()(d));
  }
  return std::make_pair(solve.eigenvectors(), sqrt_D);
}

void PT_iTEBD_step(Eigen::VectorXcd & LambdaA, Eigen::VectorXcd &LambdaB,
                   MPS_Matrix & GammaA, MPS_Matrix & GammaB, 
                   Eigen::MatrixXcd expS, const TruncatedSVD &trunc){
  //see appendix of Orus & Vidal, PRB 78, 155117 (2008)
  //Check dimensions
  if(LambdaA.size()!=GammaA.dim_d2){
    std::cerr<<"LambdaA.size()!=GammaA.dim_d2  ("<<LambdaA.size()<<" vs. "<<GammaA.dim_d2<<")!"<<std::endl;
    throw DummyException();
  }
  if(LambdaB.size()!=GammaB.dim_d2){
    std::cerr<<"LambdaB.size()!=GammaB.dim_d2  ("<<LambdaB.size()<<" vs. "<<GammaB.dim_d2<<")!"<<std::endl;
    throw DummyException();
  }
  if(LambdaA.size()!=GammaB.dim_d1){
    std::cerr<<"LambdaA.size()!=GammaB.dim_d1  ("<<LambdaA.size()<<" vs. "<<GammaB.dim_d1<<")!"<<std::endl;
    throw DummyException();
  }
  if(LambdaB.size()!=GammaA.dim_d1){
    std::cerr<<"LambdaB.size()!=GammaA.dim_d1  ("<<LambdaB.size()<<" vs. "<<GammaA.dim_d1<<")!"<<std::endl;
    throw DummyException();
  }
  int dim_i=expS.rows();
  if(expS.cols()!=dim_i){
    std::cerr<<"expS.cols()!=dim_i  ("<<expS.cols()<<" vs. "<<dim_i<<")!"<<std::endl;
    throw DummyException();
  }
  if(GammaA.dim_i!=dim_i){
    std::cerr<<"GammaA.dim_i!=dim_i  ("<<GammaA.dim_i<<" vs. "<<dim_i<<")!"<<std::endl;
    throw DummyException();
  }
  if(GammaB.dim_i!=dim_i){
    std::cerr<<"GammaB.dim_i!=dim_i  ("<<GammaB.dim_i<<" vs. "<<dim_i<<")!"<<std::endl;
    throw DummyException();
  }

std::cout<<"TEST0"<<std::endl;
  //Construct fragments:
  //Start with progagated Gamma_A LambdaA GammaB
  Eigen::MatrixXcd GLG2=Eigen::MatrixXcd::Zero(dim_i*GammaA.dim_d1, dim_i*GammaB.dim_d2);
  for(int i1=0; i1<dim_i; i1++){
    for(int i2=0; i2<dim_i; i2++){
      for(int d1=0; d1<GammaA.dim_d1; d1++){
        for(int d2=0; d2<GammaB.dim_d2; d2++){
          for(int d3=0; d3<GammaA.dim_d2; d3++){
            GLG2(i2*GammaA.dim_d1+d1, i1*GammaB.dim_d2+d2) += \
         expS(i1, i2) * GammaA(i1, d1, d3) * LambdaA(d3) * GammaB(i2,d3,d2);
          }
        }
      } 
    }
  }

std::cout<<"TEST0.5"<<std::endl;
Eigen::MatrixXcd compare=Eigen::MatrixXcd::Zero(dim_i*GammaA.dim_d1, dim_i*GammaB.dim_d2);
  for(int i1=0; i1<dim_i; i1++){
    for(int i2=0; i2<dim_i; i2++){
      for(int d1=0; d1<GammaA.dim_d1; d1++){
        for(int d2=0; d2<GammaB.dim_d2; d2++){
          compare(i1*GammaA.dim_d1+d1, i2*GammaB.dim_d2+d2) = \
LambdaB(d1)*GLG2(i1*GammaA.dim_d1+d1, i2*GammaB.dim_d2+d2)*LambdaB(d2);
        }
      } 
    }
  }
  
std::cout<<"TEST1"<<std::endl;
  //Normalize: 
  Eigen::MatrixXcd X, Xinv;
  { 
    std::pair<Eigen::MatrixXcd, Eigen::VectorXcd> Xs=PT_iTEBD_X(LambdaB, GLG2);
std::cout<<"TEST2"<<std::endl;
    X=Xs.first*Xs.second.asDiagonal();
    Xinv=Xs.second.cwiseInverse().asDiagonal()*Xs.first.adjoint();

std::cout<<"(X*Xinv-I).norm()="<<(X*Xinv-Eigen::MatrixXcd::Identity(LambdaB.size(), LambdaB.size())).norm()<<std::endl;
Eigen::MatrixXcd VR=X*X.adjoint();
Eigen::MatrixXcd R=PT_iTEBD_calc_R(LambdaB, GLG2);
Eigen::VectorXcd v(VR.rows()*VR.cols());
for(int d1=0; d1<VR.rows(); d1++){for(int d2=0; d2<VR.rows(); d2++){
v(d1*VR.cols()+d2)=VR(d1,d2);}}
Eigen::VectorXcd v2=R*v;
double norm1=std::abs(v.dot(v)); double norm2=std::abs(v2.dot(v2)); 
double scal=std::abs(v.dot(v2));    
std::cout<<"norm1="<<norm1<<" norm2="<<norm2<<" scal="<<scal<<" norm1*norm2="<<norm1*norm2<<" scal*scal="<<scal*scal<<std::endl;
  }

std::cout<<"TEST3"<<std::endl;
  Eigen::MatrixXcd YT, YTinv;
  {
    std::pair<Eigen::MatrixXcd, Eigen::VectorXcd> Ys=PT_iTEBD_Y(LambdaB, GLG2);
    std::cout<<Ys.second.transpose()<<std::endl;
//    std::cout<<Ys.first<<std::endl;
    std::cout<<"LambdaB.size()="<<LambdaB.size()<<std::endl;
std::cout<<"TEST4"<<std::endl;
    YT=(Ys.first*Ys.second.asDiagonal()).transpose();
    YTinv=(Ys.second.cwiseInverse().asDiagonal()*Ys.first.adjoint()).transpose();

std::cout<<"(YTinv*YT-I).norm()="<<(YTinv*YT-Eigen::MatrixXcd::Identity(LambdaB.size(), LambdaB.size())).norm()<<std::endl;
Eigen::MatrixXcd VL=YT.transpose()*YT.conjugate();
std::cout<<"test0"<<std::endl;
Eigen::MatrixXcd L=PT_iTEBD_calc_L(LambdaB, GLG2);
std::cout<<"test1"<<std::endl;
Eigen::VectorXcd v(VL.rows()*VL.cols());
std::cout<<"test2"<<std::endl;
for(int d1=0; d1<VL.rows(); d1++){for(int d2=0; d2<VL.rows(); d2++){
v(d1*VL.cols()+d2)=VL(d1,d2);}}
std::cout<<"test3"<<std::endl;
Eigen::VectorXcd v2=v.transpose()*L;
std::cout<<"test4"<<std::endl;
double norm1=std::abs(v.dot(v)); double norm2=std::abs(v2.dot(v2)); 
double scal=std::abs(v.dot(v2));    
std::cout<<"norm1="<<norm1<<" norm2="<<norm2<<" scal="<<scal<<" norm1*norm2="<<norm1*norm2<<" scal*scal="<<scal*scal<<std::endl;
 
  }
std::cout<<"TEST5"<<std::endl;
 
  TruncatedSVD_Return ret=trunc.compress(YT*X);

std::cout<<"TEST6"<<std::endl;
  int new_d=ret.sigma.size();
  Eigen::MatrixXcd lVXil=ret.sigma.asDiagonal()*ret.Vdagger*Xinv*LambdaB.asDiagonal();
std::cout<<"TEST7"<<std::endl;
  Eigen::MatrixXcd lYTiUl=LambdaB.asDiagonal()*YTinv*ret.U*ret.sigma.asDiagonal();
std::cout<<"TEST8"<<std::endl;

  LambdaB=ret.sigma;

  Eigen::MatrixXcd Sigma_tmp=Eigen::MatrixXcd::Zero(dim_i*new_d, dim_i*GammaB.dim_d2);
  for(int i1=0; i1<dim_i; i1++){
    for(int i2=0; i2<dim_i; i2++){
      for(int d1=0; d1<new_d; d1++){
        for(int d=0; d<GammaA.dim_d1; d++){
          for(int d2=0; d2<GammaB.dim_d2; d2++){
           Sigma_tmp(i1*new_d+d1, i2*GammaB.dim_d2+d2) += \
             lVXil(d1, d)*GLG2(i1*GammaA.dim_d1+d, i2*GammaB.dim_d2+d2);
          }
        }
      } 
    }
  }
std::cout<<"TEST9"<<std::endl;
  Eigen::MatrixXcd Sigma=Eigen::MatrixXcd::Zero(dim_i*new_d, dim_i*new_d);
  for(int i1=0; i1<dim_i; i1++){
    for(int i2=0; i2<dim_i; i2++){
      for(int d1=0; d1<new_d; d1++){
        for(int d=0; d<GammaB.dim_d2; d++){
          for(int d2=0; d2<new_d; d2++){
           Sigma(i1*new_d+d1, i2*new_d+d2) += \
             Sigma_tmp(i1*new_d+d1, i2*GammaB.dim_d2+d)*lYTiUl(d, d2);
          }
        }
      } 
    }
  }
std::cout<<"TEST10"<<std::endl;
  
  TruncatedSVD_Return ret2=trunc.compress(Sigma);
  int new_d2=ret2.sigma.size();
  LambdaA=ret2.sigma;
  GammaA.resize(dim_i, new_d, new_d2);
  for(int i=0; i<dim_i; i++){
    for(int d1=0; d1<new_d; d1++){
      for(int d2=0; d2<new_d2; d2++){
        GammaA(i,d1,d2)=ret2.U(i*new_d+d1, d2)/LambdaB(d1);
      }
    }
  }
std::cout<<"TEST11"<<std::endl;
  GammaB.resize(dim_i, new_d2, new_d);
  for(int i=0; i<dim_i; i++){
    for(int d1=0; d1<new_d2; d1++){
      for(int d2=0; d2<new_d; d2++){
        GammaB(i,d1,d2)=ret2.Vdagger(d1, i*new_d+d2)/LambdaB(d2);
      }
    }
  }

std::cout<<"TEST12"<<std::endl;
 {
  Eigen::MatrixXcd compare2=Eigen::MatrixXcd::Zero(dim_i*GammaA.dim_d1, dim_i*GammaB.dim_d2);
  for(int i1=0; i1<dim_i; i1++){
    for(int i2=0; i2<dim_i; i2++){
      for(int d1=0; d1<GammaA.dim_d1; d1++){
        for(int d2=0; d2<GammaB.dim_d2; d2++){
          for(int d3=0; d3<GammaA.dim_d2; d3++){
            compare2(i1*GammaA.dim_d1+d1, i2*GammaB.dim_d2+d2) += \
LambdaB(d1) * GammaA(i1, d1, d3) * LambdaA(d3) * GammaB(i2,d3,d2) * LambdaB(d2);
          }
        }
      } 
    }
  }
  Eigen::MatrixXcd YTinvU=YTinv*ret.U; 
  Eigen::MatrixXcd VdagXinv=ret.Vdagger*Xinv;

  Eigen::MatrixXcd tmp=Eigen::MatrixXcd::Zero(dim_i*YTinvU.rows(), dim_i*GammaB.dim_d2);
  for(int i1=0; i1<dim_i; i1++){
    for(int i2=0; i2<dim_i; i2++){
      for(int d1=0; d1<YTinvU.rows(); d1++){
        for(int d2=0; d2<GammaB.dim_d2; d2++){
          for(int d3=0; d3<GammaA.dim_d1; d3++){
            tmp(i1*GammaA.dim_d1+d1, i2*GammaB.dim_d2+d2) += \
            YTinvU(d1,d3)*compare2(i1*GammaA.dim_d1+d3, i2*GammaB.dim_d2+d2);
          }
        }
      } 
    }
  }
  Eigen::MatrixXcd tmp2=Eigen::MatrixXcd::Zero(dim_i*YTinvU.rows(), dim_i*VdagXinv.cols());
  for(int i1=0; i1<dim_i; i1++){
    for(int i2=0; i2<dim_i; i2++){
      for(int d1=0; d1<YTinvU.rows(); d1++){
        for(int d2=0; d2<VdagXinv.cols(); d2++){
          for(int d3=0; d3<GammaB.dim_d2; d3++){
            tmp2(i1*GammaA.dim_d1+d1, i2*GammaB.dim_d2+d2) += \
            tmp(i1*GammaA.dim_d1+d1, i2*GammaB.dim_d2+d3)*VdagXinv(d3,d2);
          }
        }
      } 
    }
  } 
  std::cout<<"(compare-YTinv*U*compare2*Vdagger*Xinv).norm()="<<(compare-tmp2).norm()<<std::endl;
 }
}


void PT_TEBD_step(Eigen::VectorXcd & LambdaA, const Eigen::VectorXcd &LambdaB,
                   MPS_Matrix & GammaA, MPS_Matrix & GammaB, 
                   const Eigen::MatrixXcd expS, const TruncatedSVD &trunc){
  //see appendix of Orus & Vidal, PRB 78, 155117 (2008)
  //Check dimensions
  if(LambdaA.size()!=GammaA.dim_d2){
    std::cerr<<"LambdaA.size()!=GammaA.dim_d2  ("<<LambdaA.size()<<" vs. "<<GammaA.dim_d2<<")!"<<std::endl;
    throw DummyException();
  }
  if(LambdaB.size()!=GammaB.dim_d2){
    std::cerr<<"LambdaB.size()!=GammaB.dim_d2  ("<<LambdaB.size()<<" vs. "<<GammaB.dim_d2<<")!"<<std::endl;
    throw DummyException();
  }
  if(LambdaA.size()!=GammaB.dim_d1){
    std::cerr<<"LambdaA.size()!=GammaB.dim_d1  ("<<LambdaA.size()<<" vs. "<<GammaB.dim_d1<<")!"<<std::endl;
    throw DummyException();
  }
  if(LambdaB.size()!=GammaA.dim_d1){
    std::cerr<<"LambdaB.size()!=GammaA.dim_d1  ("<<LambdaB.size()<<" vs. "<<GammaA.dim_d1<<")!"<<std::endl;
    throw DummyException();
  }
  int dim_i=expS.rows();
  if(expS.cols()!=dim_i){
    std::cerr<<"expS.cols()!=dim_i  ("<<expS.cols()<<" vs. "<<dim_i<<")!"<<std::endl;
    throw DummyException();
  }
  if(GammaA.dim_i!=dim_i){
    std::cerr<<"GammaA.dim_i!=dim_i  ("<<GammaA.dim_i<<" vs. "<<dim_i<<")!"<<std::endl;
    throw DummyException();
  }
  if(GammaB.dim_i!=dim_i){
    std::cerr<<"GammaB.dim_i!=dim_i  ("<<GammaB.dim_i<<" vs. "<<dim_i<<")!"<<std::endl;
    throw DummyException();
  }



  //Construct fragments:
  //Start with progagated Gamma_A LambdaA GammaB
  Eigen::MatrixXcd fragment=Eigen::MatrixXcd::Zero(dim_i*GammaA.dim_d1, dim_i*GammaB.dim_d2);
  for(int i1=0; i1<dim_i; i1++){
    for(int i2=0; i2<dim_i; i2++){
      for(int d1=0; d1<GammaA.dim_d1; d1++){
        for(int d2=0; d2<GammaB.dim_d2; d2++){
          for(int d3=0; d3<GammaA.dim_d2; d3++){
fragment(i2*GammaA.dim_d1+d1, i1*GammaB.dim_d2+d2) += expS(i1, i2) *  \
LambdaB(d1) * GammaA(i1, d1, d3) * LambdaA(d3) * GammaB(i2,d3,d2) * LambdaB(d2);
          }
        }
      } 
    }
  }

  TruncatedSVD_Return ret=trunc.compress(fragment);
  int new_d=ret.sigma.size();
  LambdaA=ret.sigma;

  Eigen::VectorXcd LambdaB_inv=LambdaB.cwiseInverse();
/*
  Eigen::VectorXcd LambdaB_inv(LambdaB.size());
  for(int d=0; d<LambdaB.size(); d++){
    LambdaB_inv(d)=1./LambdaB(d);
  } 
*/

  GammaA.resize(dim_i, GammaA.dim_d1, new_d);
  for(int i=0; i<dim_i; i++){
    for(int d1=0; d1<GammaA.dim_d1; d1++){
      for(int d2=0; d2<new_d; d2++){
        GammaA(i, d1, d2)=LambdaB_inv(d1)*ret.U(i*GammaA.dim_d1+d1, d2);
      }
    }
  }
  GammaB.resize(dim_i, new_d, GammaB.dim_d2); 
  for(int i=0; i<dim_i; i++){
    for(int d1=0; d1<new_d; d1++){
      for(int d2=0; d2<GammaB.dim_d2; d2++){
        GammaB(i, d1, d2)=ret.Vdagger(d1,i*GammaB.dim_d2+d2)*LambdaB_inv(d2);
      }
    }
  }
}


int main(int args, char **argv){
  Parameters param(args, argv, true);
  
  std::string write_PT=param.get_as_string_check("write_PT");
  double dict_zero=param.get_as_double("dict_zero",0);
  TruncatedSVD trunc(param);

  DiagBB diagBB;
  {std::string Gaussian_prefix=param.get_as_string("Gaussian_prefix", "Boson");
   diagBB.setup(param, Gaussian_prefix);}


  if(!diagBB.is_set_up()){
    std::cerr<<"DiagBB not set up!"<<std::endl;
    throw DummyException();
  }
  int N=diagBB.get_dim();  
  int NL=N*N;

  TimeGrid tgrid(param);
  int n_mem=tgrid.n_mem;
  if(tgrid.n_mem<2){
    std::cerr<<"tgrid.n_mem<2!"<<std::endl;
    throw DummyException();
  }

  Eigen::VectorXcd LambdaA(1); LambdaA(0)=1;
  Eigen::VectorXcd LambdaB(1); LambdaB(0)=1;
  MPS_Matrix GammaA(NL+1, 1,1); GammaA.fill(1.);
  MPS_Matrix GammaB(NL+1, 1,1); GammaB.fill(1.);

  for(int n=n_mem-1; n>0; n--){ 
    std::cout<<"n="<<n<<"/"<<n_mem<<" dim_A="<<LambdaA.size()<<" dim_B="<<LambdaB.size()<<std::endl;
//std::cout<<"TEST0 n="<<n<<std::endl;
    Eigen::MatrixXcd expS=Eigen::MatrixXcd::Ones(NL+1,NL+1);
    expS.block(0,0,NL,NL)=diagBB.calculate_expS(n,tgrid.dt);
    if(n%2==0){
//      PT_iTEBD_step(LambdaA, LambdaB, GammaA, GammaB, expS, trunc);
      PT_TEBD_step(LambdaA, LambdaB, GammaA, GammaB, expS, trunc);
    }else{
//      PT_iTEBD_step(LambdaB, LambdaA, GammaB, GammaA, expS, trunc);
      PT_TEBD_step(LambdaB, LambdaA, GammaB, GammaA, expS, trunc);
    }
  }

  //First time step;
  std::cout<<"n="<<0<<"/"<<n_mem<<" dim_A="<<LambdaA.size()<<" dim_B="<<LambdaB.size()<<std::endl;
  MPS_Matrix f(NL+1, GammaA.dim_d1, GammaB.dim_d2); 
  f.set_zero();
  Eigen::MatrixXcd expS=Eigen::MatrixXcd::Ones(NL+1,NL+1);
  expS.block(0,0,NL,NL)=diagBB.calculate_expS(0,tgrid.dt);
  for(int i=0; i<NL+1; i++){
    f.get_Matrix_d1_d2(i) = expS(i,i) * \
                           GammaA.get_Matrix_d1_d2(i)*LambdaA.asDiagonal()*\
                           GammaB.get_Matrix_d1_d2(i)*LambdaB.asDiagonal();
  }
/*
  std::cout<<"GammaA.dim_d1="<<GammaA.dim_d1;
  std::cout<<" GammaA.dim_d2="<<GammaA.dim_d2;
  std::cout<<" GammaB.dim_d1="<<GammaB.dim_d1;
  std::cout<<" GammaB.dim_d2="<<GammaB.dim_d2<<std::endl;
*/

  Eigen::MatrixXcd frag0 = f.get_Matrix_d1_d2(NL);

  Eigen::VectorXcd vr, vl;
  { Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(frag0);
//    std::cout<<"Right eigenvalues: "<<solver.eigenvalues().transpose()<<std::endl;
    int dim_d=LambdaB.size();
    std::vector<std::pair<int, std::complex<double> > > srt(dim_d);
    for(int d=0; d<dim_d; d++){
      srt[d].first=d; srt[d].second=solver.eigenvalues()(d);
    }
    std::sort(srt.begin(), srt.end(), [] (const std::pair<int, std::complex<double> > & a, const std::pair<int, std::complex<double> > &b){ return std::abs(a.second)>std::abs(b.second);});
    std::cout<<"Largest eigenvalue: "<<srt[0].second<<" (nr="<<srt[0].first<<")"<<std::endl;
    vr=solver.eigenvectors().col(srt[0].first);
  }

  { Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(frag0.transpose());
//    std::cout<<"Left eigenvalues: "<<solver.eigenvalues().transpose()<<std::endl;
    int dim_d=LambdaB.size();
    std::vector<std::pair<int, std::complex<double> > > srt(dim_d);
    for(int d=0; d<dim_d; d++){
      srt[d].first=d; srt[d].second=solver.eigenvalues()(d);
    }
    std::sort(srt.begin(), srt.end(), [] (const std::pair<int, std::complex<double> > & a, const std::pair<int, std::complex<double> > &b){ return std::abs(a.second)>std::abs(b.second);});
    std::cout<<"Largest eigenvalue: "<<srt[0].second<<" (nr="<<srt[0].first<<")"<<std::endl;
    vl=solver.eigenvectors().col(srt[0].first).transpose();
  }

  int buffer_blocksize=param.get_as_int("buffer_blocksize",-1);

  if(param.get_as_bool("write_as_PTB",false)){
    ProcessTensorBuffer PTB(ProcessTensorBufferSpec(write_PT, buffer_blocksize));
    PTB.resize(tgrid.n_tot);
    for(int n=0; n<tgrid.n_tot; n++){
      ProcessTensorElement &e = PTB.get(n);
      e.clear();
      e.accessor.dict.set_default_diag(N);
      e.M.resize(NL,f.dim_d1,f.dim_d2);
      for(int i=0; i<NL; i++){
        e.M.get_Matrix_d1_d2(i)=f.get_Matrix_d1_d2(i);
      }
    } 
    std::complex<double> prod=vl.transpose()*vr;
    vr/=prod;
    std::cout<<"vl.transpose()*vr="<<vl.transpose()*vr<<std::endl;
    std::cout<<"vl.transpose()*frag0*vr="<<vl.transpose()*frag0*vr<<std::endl;
    PTB.get(0).M.inner_multiply_left(vl.transpose());
    PTB.get(tgrid.n_tot-1).M.inner_multiply_right(vr);
    PTB.expand_DiagBB(diagBB, dict_zero);
    PTB.calculate_closures();

  }else{ //write as ProcessTensorRepeat
    ProcessTensorBuffer PTB;
    PTB.resize(3);
    for(int n=0; n<3; n++){
      ProcessTensorElement &e = PTB.get(n);
      e.clear();
      e.accessor.dict.set_default_diag(N);
      e.M.resize(NL,f.dim_d1,f.dim_d2);
      for(int i=0; i<NL; i++){
        e.M.get_Matrix_d1_d2(i)=f.get_Matrix_d1_d2(i);
      }
    } 
    std::complex<double> prod=vl.transpose()*vr;
    vr/=prod;
    std::cout<<"vl.transpose()*vr="<<vl.transpose()*vr<<std::endl;
    std::cout<<"vl.transpose()*frag0*vr="<<vl.transpose()*frag0*vr<<std::endl;
    PTB.get(0).M.inner_multiply_left(vl.transpose());
    PTB.get(2).M.inner_multiply_right(vr);
    PTB.expand_DiagBB(diagBB, dict_zero);
    PTB.calculate_closures();
  
    ProcessTensorRepeat PTR;
    PTR.set_specs(write_PT, buffer_blocksize);
    PTR.initial.resize(1);
    PTR.initial.get(0)=PTB.get(0);
    PTR.repeated.resize(1);
    PTR.repeated.get(0)=PTB.get(1);
  }
 
  return 0;
}
