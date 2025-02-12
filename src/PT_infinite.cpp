#include "ACE.hpp"
#include "PT_infinite.hpp"
#include "ProcessTensorBuffer.hpp"
#include "ProcessTensorRepeat.hpp"
#include "Largest_EV.hpp"
#include "Timings.hpp"

namespace ACE{


Eigen::MatrixXcd PT_iTEBD_calc_R(const Eigen::VectorXcd & LambdaB, const Eigen::MatrixXcd &GLG2){
  int dim_d=LambdaB.size();
  int dim_i=GLG2.rows()/dim_d;
  Eigen::MatrixXcd R=Eigen::MatrixXcd::Zero(dim_d*dim_d, dim_d*dim_d);
  for(int d3=0; d3<dim_d; d3++){
    for(int d4=0; d4<dim_d; d4++){
      for(int j=0; j<dim_i; j++){
        for(int d1=0; d1<dim_d; d1++){
          for(int d2=0; d2<dim_d; d2++){
            for(int i=0; i<dim_i; i++){
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
  for(int d3=0; d3<dim_d; d3++){
    for(int d4=0; d4<dim_d; d4++){
      for(int j=0; j<dim_i; j++){
        for(int d1=0; d1<dim_d; d1++){
          for(int d2=0; d2<dim_d; d2++){
            for(int i=0; i<dim_i; i++){
              L(d1*dim_d+d2, d3*dim_d+d4)+=  \
LambdaB(d1)*GLG2(i*dim_d+d1, j*dim_d+d3)*std::conj(LambdaB(d2)*GLG2(i*dim_d+d2,j*dim_d+d4));
            }
          }
        }
      }
    }
  }
  return L;
}
Eigen::VectorXcd PT_iTEBD_R_times(const Eigen::VectorXcd & LambdaB, const Eigen::MatrixXcd &GLG2, const Eigen::VectorXcd &vr){

  int dim_d=LambdaB.size();
  int dim_i=GLG2.rows()/dim_d;

/*
  typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXcdRM;
  MatrixXcdRM VR2=MatrixXcdRM::Zero(dim_d, dim_d);
  MatrixXcdRM VR_ = LambdaB.asDiagonal() * \
         Eigen::Map<const MatrixXcdRM>(&vr(0), dim_d, dim_d) * \
                         LambdaB.conjugate().asDiagonal();
//  Eigen::MatrixXcd VR2=Eigen::MatrixXcd::Zero(dim_d, dim_d);
//  Eigen::MatrixXcd VR_ = LambdaB.asDiagonal() * \
                         Eigen::Map<const Eigen::MatrixXcd>(&vr(0), dim_d, dim_d) * \
                         LambdaB.conjugate().asDiagonal();

  for(int i=0; i<dim_i; i++){
   for(int j=0; j<dim_i; j++){
      VR2.noalias()+= 
          Eigen::Map<const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>, 0, Eigen::OuterStride<> > (&GLG2(i*dim_d,j*dim_d), dim_d, dim_d, Eigen::OuterStride<>(dim_i*dim_d)) * VR_ * \
          Eigen::Map<const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::OuterStride<> > (&GLG2(i*dim_d,j*dim_d), dim_d, dim_d, Eigen::OuterStride<>(dim_i*dim_d)).conjugate();

    }
  }
  return Eigen::Map<Eigen::VectorXcd>(&VR_(0,0), dim_d*dim_d);
*/


  //Note: Matrix is transposed as matrix-vector mapping is ColMajor:
  Eigen::MatrixXcd VR2=Eigen::MatrixXcd::Zero(dim_d, dim_d);
  Eigen::MatrixXcd VR_ = LambdaB.conjugate().asDiagonal() * \
                         Eigen::Map<const Eigen::MatrixXcd>(&vr(0), dim_d, dim_d) * \
                         LambdaB.asDiagonal();

  for(int i=0; i<dim_i; i++){
   for(int j=0; j<dim_i; j++){
      VR2.noalias()+= 
        Eigen::Map<const Eigen::MatrixXcd, 0, Eigen::OuterStride<> > (&GLG2(i*dim_d,j*dim_d), dim_d, dim_d, Eigen::OuterStride<>(dim_i*dim_d)).conjugate()  * \
                                 VR_ * \
        Eigen::Map<const Eigen::MatrixXcd, 0, Eigen::OuterStride<> > (&GLG2(i*dim_d,j*dim_d), dim_d, dim_d, Eigen::OuterStride<>(dim_i*dim_d)).transpose();
//          Eigen::Map<const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::OuterStride<> > (&GLG2(i*dim_d,j*dim_d), dim_d, dim_d, Eigen::OuterStride<>(dim_i*dim_d)).adjoint()  * \
                                 VR_ * \
          Eigen::Map<const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>, 0, Eigen::OuterStride<> > (&GLG2(i*dim_d,j*dim_d), dim_d, dim_d, Eigen::OuterStride<>(dim_i*dim_d)).transpose();
    }
  }
  return Eigen::Map<Eigen::VectorXcd>(&VR2(0,0), dim_d*dim_d);
/* 
  Eigen::VectorXcd vr2 = Eigen::VectorXcd::Zero(dim_d*dim_d);
  for(int i=0; i<dim_i; i++){
   for(int j=0; j<dim_i; j++){
      Eigen::VectorXcd tmp=Eigen::VectorXcd::Zero(dim_d*dim_d);
      for(int d1=0; d1<dim_d; d1++){
        for(int d2=0; d2<dim_d; d2++){
          for(int d=0; d<dim_d; d++){
            tmp(d1*dim_d+d2)+=GLG2(i*dim_d+d1, j*dim_d+d)*LambdaB(d)*vr(d*dim_d+d2);
          }
        }
      }
      for(int d1=0; d1<dim_d; d1++){
        for(int d2=0; d2<dim_d; d2++){
          for(int d=0; d<dim_d; d++){
            vr2(d1*dim_d+d2)+=std::conj(GLG2(i*dim_d+d2,j*dim_d+d)*LambdaB(d))*tmp(d1*dim_d+d);
          }
        }
      }
    }
  }
  return vr2;
*/
}

Eigen::VectorXcd PT_iTEBD_LT_times(const Eigen::VectorXcd & LambdaB, const Eigen::MatrixXcd &GLG2, const Eigen::VectorXcd &vl){
  int dim_d=LambdaB.size();
  int dim_i=GLG2.rows()/dim_d;

/*
  Eigen::MatrixXcd VL2=Eigen::MatrixXcd::Zero(dim_d, dim_d);
  Eigen::MatrixXcd VL_ = LambdaB.asDiagonal() * \
                         Eigen::Map<const Eigen::MatrixXcd>(&vl(0), dim_d, dim_d) * \
                         LambdaB.conjugate().asDiagonal();
  for(int i=0; i<dim_i; i++){
    for(int j=0; j<dim_i; j++){
      VL2.noalias()+= 
          Eigen::Map<const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::OuterStride<> > (&GLG2(i*dim_d,j*dim_d), dim_d, dim_d, Eigen::OuterStride<>(dim_i*dim_d)) * VL_ * \
          Eigen::Map<const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>, 0, Eigen::OuterStride<> > (&GLG2(i*dim_d,j*dim_d), dim_d, dim_d, Eigen::OuterStride<>(dim_i*dim_d)).conjugate();
    }
  }
  return Eigen::Map<Eigen::VectorXcd>(&VL_(0,0), dim_d*dim_d);
*/

//Note: Matrices VL* are transposed (ColumnMajor)
  Eigen::MatrixXcd VL2=Eigen::MatrixXcd::Zero(dim_d, dim_d);
  Eigen::MatrixXcd VL_ = LambdaB.conjugate().asDiagonal() * \
                         Eigen::Map<const Eigen::MatrixXcd>(&vl(0), dim_d, dim_d) * \
                         LambdaB.asDiagonal();
  for(int i=0; i<dim_i; i++){
    for(int j=0; j<dim_i; j++){
      VL2.noalias()+= 
          Eigen::Map<const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>, 0, Eigen::OuterStride<> > (&GLG2(i*dim_d,j*dim_d), dim_d, dim_d, Eigen::OuterStride<>(dim_i*dim_d)).adjoint() *\
                                        VL_ * \
          Eigen::Map<const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::OuterStride<> > (&GLG2(i*dim_d,j*dim_d), dim_d, dim_d, Eigen::OuterStride<>(dim_i*dim_d)).transpose();
    }
  }
  return Eigen::Map<Eigen::VectorXcd>(&VL2(0,0), dim_d*dim_d);
/*
  Eigen::VectorXcd vl2=Eigen::VectorXcd::Zero(dim_d*dim_d);
  for(int i=0; i<dim_i; i++){
    for(int j=0; j<dim_i; j++){
      Eigen::VectorXcd tmp=Eigen::VectorXcd::Zero(dim_d*dim_d);
      for(int d1=0; d1<dim_d; d1++){
        for(int d2=0; d2<dim_d; d2++){
          for(int d=0; d<dim_d; d++){
            tmp(d1*dim_d+d2)+=vl(d*dim_d+d2)*LambdaB(d)*GLG2(i*dim_d+d, j*dim_d+d1);
          }
        }
      }
      for(int d1=0; d1<dim_d; d1++){
        for(int d2=0; d2<dim_d; d2++){
          for(int d=0; d<dim_d; d++){
            vl2(d1*dim_d+d2)+=tmp(d1*dim_d+d)*std::conj(LambdaB(d)*GLG2(i*dim_d+d,j*dim_d+d2));
          }
        }
      }
    }
  }
  return vl2;
*/
}
 

std::pair<Eigen::MatrixXcd, Eigen::MatrixXcd> PT_iTEBD_X(const Eigen::VectorXcd & LambdaB, const Eigen::MatrixXcd &GLG2, const infinite_normalize_specs & specs){

  int dim_d=LambdaB.size();
  int dim_i=GLG2.rows()/dim_d;

   //Largest eigenvalue -> matrix VR:
  Eigen::MatrixXcd VR=Eigen::MatrixXcd::Identity(dim_d, dim_d)/sqrt((double)dim_d);

/*
  //start with good guess?
  {  
    Eigen::MatrixXcd R=Eigen::Map<const Eigen::MatrixXcd>(&GLG2((dim_i-1)*dim_d,(dim_i-1)*dim_d), dim_d, dim_d)*LambdaB.asDiagonal();
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solve_initial(R);
    Eigen::MatrixXcd U=solve_initial.eigenvectors();  
//    VR=U*U.adjoint();
    int imax=0; for(int i=1; i<dim_d; i++){if(std::real(solve_initial.eigenvalues()(i))>std::real(solve_initial.eigenvalues()(imax))){imax=i;}}
    VR=solve_initial.eigenvectors().col(imax)*(solve_initial.eigenvectors().col(imax)).adjoint();

    VR/=VR.norm();

    Eigen::VectorXcd vr(dim_d*dim_d);
    for(int d1=0; d1<dim_d; d1++){
      for(int d2=0; d2<dim_d; d2++){
        vr(d1*dim_d+d2)=VR(d1,d2);
      }
    }
    //Check how close to an EV we are:
    Eigen::VectorXcd guess2=PT_iTEBD_R_times(LambdaB, GLG2, vr);
    std::complex<double> guess_max_eval=vr.dot(guess2);
    std::cout<<"|guess2/guess_max_eval-guess|="<<(guess2/guess_max_eval-vr).norm()<<std::endl;
    std::cout<<"guess_max_eval="<<guess_max_eval<<std::endl;
std::cout<<"guess:"<<vr.transpose()<<std::endl;
  }
*/

  std::complex<double> max_eval=0;
  if(!specs.use_iter){ //Brute-force solution of eigenvalue problem
//    time_point beforeR=now();
    Eigen::MatrixXcd R=PT_iTEBD_calc_R(LambdaB, GLG2);
//    std::cout<<"Time to calculate R: "<<time_diff(now()-beforeR)<<std::endl;

    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solveR(R);
    int imax=0;
    for(int i=1; i<dim_d*dim_d; i++){
      if(std::real(solveR.eigenvalues()(i))>std::real(solveR.eigenvalues()(imax))){
        imax=i;
      }
    }
    for(int d1=0; d1<dim_d; d1++){ 
      for(int d2=0; d2<dim_d; d2++){
        VR(d1,d2)=solveR.eigenvectors()(d1*dim_d+d2,imax);
      }
    } 
    max_eval=solveR.eigenvalues()(imax);

  }else if(!specs.use_Arnoldi){ //Power iteration
    Eigen::VectorXcd vr=Eigen::Map<Eigen::VectorXcd>(&VR(0,0), dim_d*dim_d);
    for(int it=1; it<=specs.iter; it++){
      Eigen::VectorXcd vr2=PT_iTEBD_R_times(LambdaB, GLG2, vr);
      max_eval=vr.dot(vr2);
      double diff=(vr2/max_eval-vr).norm();
      vr.swap(vr2); vr.normalize();
      if(it>1 && diff<specs.eps){
        std::cout<<"R: converged to "<<diff<<"<"<<specs.eps<<" at iteration "<<it<<"/"<<specs.iter<<std::endl;
        break;
      } 
/*
      double nrm=VR2.norm();
      VR2/=nrm;
      double diff=(VR-VR2).norm();
      if(it%100==0){
        std::cout<<"R: iter="<<it<<" diff="<<diff<<std::endl;
      }
      VR.swap(VR2);
      max_eval=nrm;
      if(it>1 && diff<epsilon){
        std::cout<<"R: converged to "<<diff<<"<"<<epsilon<<" at iteration "<<it<<"/"<<iter<<std::endl;
        break;
      }
*/
    }
    VR=Eigen::Map<const Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic, Eigen::RowMajor> >(&vr(0), dim_d, dim_d);

  }else{ //Arnoldi
    Eigen::VectorXcd vr(dim_d*dim_d);
    for(int d1=0; d1<dim_d; d1++){
      for(int d2=0; d2<dim_d; d2++){
        vr(d1*dim_d+d2)=VR(d1,d2);
      }
    }
//std::cout<<"LmabdaB:"<<LambdaB.transpose()<<std::endl;
//std::cout<<"GLG2:"<<std::endl<<GLG2<<std::endl;
    max_eval=Largest_EV_Arnoldi_BLAS(vr, specs.iter, specs.eps, 
       [&LambdaB,&GLG2](const Eigen::VectorXcd &v){
         return PT_iTEBD_R_times(LambdaB, GLG2, v);
       }, 1);
    for(int d1=0; d1<dim_d; d1++){
      for(int d2=0; d2<dim_d; d2++){
        VR(d1,d2)=vr(d1*dim_d+d2);
      }
    }

    //Check how close to an EV we are:
    Eigen::VectorXcd vr2=PT_iTEBD_R_times(LambdaB, GLG2, vr);
    std::cout<<"|VR2/max_eval-VR|="<<(vr2/max_eval-vr).norm()<<std::endl;

std::cout<<"max_eval="<<max_eval<<std::endl;
//std::cout<<"vr:"<<vr.transpose()<<std::endl;
/*
std::cout<<"max_eval="<<max_eval<<std::endl;
    //Compare:
    Eigen::MatrixXcd R=PT_iTEBD_calc_R(LambdaB, GLG2);
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solveR(R);
    int imax=0; for(int i=1; i<dim_d*dim_d; i++){if(std::real(solveR.eigenvalues()(i))>std::real(solveR.eigenvalues()(imax))){imax=i;}}
    Eigen::VectorXcd compare=solveR.eigenvectors().col(imax);
    std::cout<<"1-|compare*vr|="<<1.-std::abs(compare.dot(vr))<<std::endl;
//std::cout<<"compare="<<compare.transpose()<<std::endl;
//std::cout<<"max_eval2="<<solveR.eigenvalues()(imax)<<std::endl;
//std::cout<<"all eigenvalues="<<solveR.eigenvalues().transpose()<<std::endl;

    Eigen::MatrixXcd R2(dim_d*dim_d, dim_d*dim_d);
    for(int D=0; D<dim_d*dim_d; D++){
      Eigen::VectorXcd tmp=Eigen::VectorXcd::Zero(dim_d*dim_d); tmp(D)=1.;
      R2.col(D)=PT_iTEBD_R_times(LambdaB, GLG2, tmp);
    }
std::cout<<"|R2-R|="<<(R2-R).norm()<<std::endl;
*/
  }
/*
  std::cout<<"|VR-VR.adjoint()|="<<(VR-VR.adjoint()).norm()<<std::endl;
std::cout<<"R:"<<std::endl<<R<<std::endl;
std::cout<<"VR:"<<std::endl<<VR<<std::endl;
std::cout<<"R eigenvalues"<<solveR.eigenvalues().transpose()<<std::endl;
std::cout<<"max eval="<<max_eval<<" at imax="<<imax<<std::endl;
*/
  //Decompose VR:
  std::cout<<"|VR-VR.adjoint()|="<<(VR-VR.adjoint()).norm()<<std::endl;
  Eigen::MatrixXcd input=0.5*(VR+VR.adjoint()); //VR+1e-8*Eigen::MatrixXcd::Identity(dim_d,dim_d);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solve(input);
  Eigen::VectorXcd sqrt_D(dim_d);

  int max_abs_i=0;
  for(int d=1; d<dim_d; d++){ if(std::abs(solve.eigenvalues()(d))>std::abs(solve.eigenvalues()(max_abs_i))){ max_abs_i=d; } }
  
  for(int d=0; d<dim_d; d++){  
    std::complex<double> eval=solve.eigenvalues()(d)/solve.eigenvalues()(max_abs_i);
    if(std::abs(eval)<1e-15){
      std::cerr<<"PT_iTEBD_X: not invertible: "<<solve.eigenvalues().transpose()<<std::endl;
      throw DummyException();
    }
    sqrt_D(d)=sqrt(eval);
  }
  std::cout<<"X: reconstruct VR: "<<(solve.eigenvectors()*sqrt_D.asDiagonal()*sqrt_D.asDiagonal()*solve.eigenvectors().adjoint()-input).norm()<<std::endl;
  return {solve.eigenvectors()*sqrt_D.asDiagonal(),
          sqrt_D.cwiseInverse().asDiagonal()*solve.eigenvectors().adjoint()};

/*
  Eigen::MatrixXcd input=VR+1e-8*Eigen::MatrixXcd::Ones(dim_d,dim_d);
  Eigen::LDLT<Eigen::MatrixXcd> llt(input);
  Eigen::MatrixXcd l=llt.matrixL();
  Eigen::MatrixXcd pl=llt.transpositionsP().inverse()*l;
  Eigen::MatrixXcd lp=l.triangularView<Eigen::Lower>().solve(Eigen::MatrixXcd::Identity(dim_d,dim_d));
  lp=lp*llt.transpositionsP().inverse();
  Eigen::VectorXcd Dsqrt=llt.vectorD().cwiseSqrt();
  
  std::cout<<"LDLT: "<<((pl*Dsqrt.asDiagonal())*(Dsqrt.asDiagonal()*pl.adjoint())-input).norm()<<std::endl;
  std::cout<<"LDLT_rec: "<<(llt.reconstructedMatrix()-input).norm()<<std::endl;
  std::cout<<"|lp*pl-1|="<<(lp*pl-Eigen::MatrixXcd::Identity(dim_d,dim_d)).norm()<<std::endl;
  std::cout<<"|(sD^{-1}lp)(pl*sD)-1|="<<(Dsqrt.cwiseInverse().asDiagonal()*lp*pl*Dsqrt.asDiagonal()-Eigen::MatrixXcd::Identity(dim_d,dim_d)).norm()<<std::endl;
  return {pl*Dsqrt.asDiagonal(), Dsqrt.cwiseInverse().asDiagonal()*lp};
*/
}

std::pair<Eigen::MatrixXcd, Eigen::MatrixXcd> PT_iTEBD_Y(const Eigen::VectorXcd & LambdaB, const Eigen::MatrixXcd &GLG2, const infinite_normalize_specs & specs){

  int dim_d=LambdaB.size();
  int dim_i=GLG2.rows()/dim_d;

   //Largest eigenvalue -> matrix VL:
  Eigen::MatrixXcd VL=Eigen::MatrixXcd::Identity(dim_d, dim_d)/sqrt((double)dim_d);

  std::complex<double> max_eval=0;
  if(!specs.use_iter){ //Brute-force solution of eigenvalue problem
    Eigen::MatrixXcd L=PT_iTEBD_calc_L(LambdaB, GLG2);
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solveL(L.transpose());
    int imax=0;
    for(int i=1; i<dim_d*dim_d; i++){
      if(std::real(solveL.eigenvalues()(i))>std::real(solveL.eigenvalues()(imax))){
        imax=i;
      }
    }
    for(int d1=0; d1<dim_d; d1++){ 
      for(int d2=0; d2<dim_d; d2++){
        VL(d1,d2)=solveL.eigenvectors()(d1*dim_d+d2,imax);
      }
    } 
    max_eval=solveL.eigenvalues()(imax);

  }else if(!specs.use_Arnoldi){ //Power iteration
    Eigen::VectorXcd vl=Eigen::Map<Eigen::VectorXcd>(&VL(0,0), dim_d*dim_d);
    for(int it=1; it<=specs.iter; it++){
      Eigen::VectorXcd vl2=PT_iTEBD_LT_times(LambdaB, GLG2, vl);
      max_eval=vl.dot(vl2);
      double diff=(vl2/max_eval-vl).norm();
      vl.swap(vl2); vl.normalize();
      if(it>1 && diff<specs.eps){
        std::cout<<"L: converged to "<<diff<<"<"<<specs.eps<<" at iteration "<<it<<"/"<<specs.iter<<std::endl;
        break;
      } 
/*
      Eigen::MatrixXcd VL2=PT_iTEBD_LT_times(LambdaB, GLG2, VL);

      double nrm=VL2.norm();
      VL2/=nrm;
      double diff=(VL-VL2).norm();
      if(it%100==0){
        std::cout<<"L: iter="<<it<<" diff="<<diff<<std::endl;
      }
      VL.swap(VL2);
      max_eval=nrm;
      if(it>1 && diff<epsilon){
        std::cout<<"L: converged to "<<diff<<"<"<<epsilon<<" at iteration "<<it<<"/"<<iter<<std::endl;
        break;
      }
*/
    }
    VL=Eigen::Map<const Eigen::MatrixXcd>(&vl(0), dim_d, dim_d);

  }else{ //Arnoldi
    Eigen::VectorXcd vl(dim_d*dim_d);
    for(int d1=0; d1<dim_d; d1++){
      for(int d2=0; d2<dim_d; d2++){
        vl(d1*dim_d+d2)=VL(d1,d2);
      }
    }
    max_eval=Largest_EV_Arnoldi_BLAS(vl, specs.iter, specs.eps, 
       [&LambdaB,&GLG2](const Eigen::VectorXcd &v){
         return PT_iTEBD_LT_times(LambdaB, GLG2, v);
       }, 1);
    for(int d1=0; d1<dim_d; d1++){
      for(int d2=0; d2<dim_d; d2++){
        VL(d1,d2)=vl(d1*dim_d+d2);
      }
    }

    //Check how close to an EV we are:
    Eigen::VectorXcd vl2=PT_iTEBD_LT_times(LambdaB, GLG2, vl);
    std::cout<<"|VL2/max_eval-VL|="<<(vl2/max_eval-vl).norm()<<std::endl;

std::cout<<"max_eval="<<max_eval<<std::endl;
/*
std::cout<<"max_eval="<<max_eval<<std::endl;
    //Compare:
    Eigen::MatrixXcd L=PT_iTEBD_calc_L(LambdaB, GLG2);
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solveL(L.transpose());
    int imax=0; for(int i=1; i<dim_d*dim_d; i++){if(std::real(solveL.eigenvalues()(i))>std::real(solveL.eigenvalues()(imax))){imax=i;}}
    Eigen::VectorXcd compare=solveL.eigenvectors().col(imax);
    std::cout<<"1-|compare*vl|="<<1.-std::abs(compare.dot(vl))<<std::endl;
//std::cout<<"compare="<<compare.transpose()<<std::endl;
//std::cout<<"max_eval2="<<solveL.eigenvalues()(imax)<<std::endl;
//std::cout<<"all eigenvalues="<<solveL.eigenvalues().transpose()<<std::endl;

    Eigen::MatrixXcd L2(dim_d*dim_d, dim_d*dim_d);
    for(int D=0; D<dim_d*dim_d; D++){
      Eigen::VectorXcd tmp=Eigen::VectorXcd::Zero(dim_d*dim_d); tmp(D)=1.;
      L2.col(D)=PT_iTEBD_LT_times(LambdaB, GLG2, tmp);
    }
std::cout<<"|L2-L|="<<(L2-L).norm()<<std::endl;
std::cout<<"|L2-LT|="<<(L2-L.transpose()).norm()<<std::endl;
*/
  }


  //Decompose VL:
  std::cout<<"|VL-VL.adjoint()|="<<(VL-VL.adjoint()).norm()<<std::endl;
  Eigen::MatrixXcd input=0.5*(VL+VL.adjoint()); //VL+1e-8*Eigen::MatrixXcd::Identity(dim_d,dim_d);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solve(input);
  Eigen::VectorXcd sqrt_D(dim_d);

  int max_abs_i=0;
  for(int d=1; d<dim_d; d++){ if(std::abs(solve.eigenvalues()(d))>std::abs(solve.eigenvalues()(max_abs_i))){ max_abs_i=d; } }
  
  for(int d=0; d<dim_d; d++){  
    std::complex<double> eval=solve.eigenvalues()(d)/solve.eigenvalues()(max_abs_i);
    if(std::abs(eval)<1e-15){
      std::cerr<<"PT_iTEBD_Y: not invertible: "<<solve.eigenvalues().transpose()<<std::endl;
      throw DummyException();
    }
    sqrt_D(d)=sqrt(eval);
  }

  std::cout<<"Y: reconstruct VL: "<<(solve.eigenvectors()*sqrt_D.asDiagonal()*sqrt_D.asDiagonal()*solve.eigenvectors().adjoint()-input).norm()<<std::endl;
//  return {solve.eigenvectors(), sqrt_D, max_eval};
  return {solve.eigenvectors()*sqrt_D.asDiagonal(),
          sqrt_D.cwiseInverse().asDiagonal()*solve.eigenvectors().adjoint()};
}

void PT_iTEBD_step(Eigen::VectorXcd & LambdaA, Eigen::VectorXcd &LambdaB,
                   MPS_Matrix & GammaA, MPS_Matrix & GammaB, 
                   Eigen::MatrixXcd expS, const TruncatedSVD &trunc,
                   const infinite_normalize_specs & specs){
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
  Eigen::MatrixXcd GLG2=Eigen::MatrixXcd::Zero(dim_i*GammaA.dim_d1, dim_i*GammaB.dim_d2);
  for(int i1=0; i1<dim_i; i1++){
    for(int i2=0; i2<dim_i; i2++){
      Eigen::Map<Eigen::MatrixXcd, 0, Eigen::OuterStride<> > (&GLG2(i2*GammaA.dim_d1,i1*GammaB.dim_d2), GammaA.dim_d1, GammaB.dim_d2, Eigen::OuterStride<>(dim_i*GammaA.dim_d1)) = expS(i1, i2) * \
GammaA.get_Matrix_d1_d2(i1) * LambdaA.asDiagonal() * GammaB.get_Matrix_d1_d2(i2);
/*
      for(int d1=0; d1<GammaA.dim_d1; d1++){
        for(int d2=0; d2<GammaB.dim_d2; d2++){
          for(int d3=0; d3<GammaA.dim_d2; d3++){
            GLG2(i2*GammaA.dim_d1+d1, i1*GammaB.dim_d2+d2) += \
         expS(i1, i2) * GammaA(i1, d1, d3) * LambdaA(d3) * GammaB(i2,d3,d2);
          }
        }
      } 
*/
    }
  }

// for(int normalize_loop=0; normalize_loop<2; normalize_loop++)
 {
  //Normalize: 
  Eigen::MatrixXcd X, Xinv;
  { 
    time_point norm_start=now();
    std::tie(X, Xinv)=PT_iTEBD_X(LambdaB, GLG2, specs);
    std::cout<<"Right canonicalize: runtime="<<time_diff(now()-norm_start)<<"ms"<<std::endl;
  }

  Eigen::MatrixXcd YT, YTinv;
  {
    
    time_point norm_start=now();
    Eigen::MatrixXcd Y,Yinv;
    std::tie(Y, Yinv)=PT_iTEBD_Y(LambdaB, GLG2, specs);
    YT=Y.transpose(); YTinv=Yinv.transpose();
    std::cout<<"Left canonicalize: runtime="<<time_diff(now()-norm_start)<<"ms"<<std::endl;
  }
 
  time_point SVD_start=now();
  TruncatedSVD_Return ret=trunc.compress(YT*LambdaB.asDiagonal()*X);
  std::cout<<"first SVD: runtime="<<time_diff(now()-SVD_start)<<"ms"<<std::endl;


  //Build new GLG2 matrix
  int new_d=ret.sigma.size();
  LambdaB=ret.sigma;
  Eigen::MatrixXcd VXi=ret.Vdagger*Xinv;
  Eigen::MatrixXcd YTiU=YTinv*ret.U;

  Eigen::MatrixXcd GLG2_tmp=Eigen::MatrixXcd::Zero(dim_i*new_d, dim_i*GammaB.dim_d2);
  for(int i1=0; i1<dim_i; i1++){
    for(int i2=0; i2<dim_i; i2++){
      for(int d1=0; d1<new_d; d1++){
        for(int d=0; d<GammaA.dim_d1; d++){
          for(int d2=0; d2<GammaB.dim_d2; d2++){
           GLG2_tmp(i1*new_d+d1, i2*GammaB.dim_d2+d2) += \
             VXi(d1, d)*GLG2(i1*GammaA.dim_d1+d, i2*GammaB.dim_d2+d2);
          }
        }
      } 
    }
  }
  GLG2=Eigen::MatrixXcd::Zero(dim_i*new_d, dim_i*new_d);
  for(int i1=0; i1<dim_i; i1++){
    for(int i2=0; i2<dim_i; i2++){
      for(int d1=0; d1<new_d; d1++){
        for(int d=0; d<GammaB.dim_d2; d++){
          for(int d2=0; d2<new_d; d2++){
           GLG2(i1*new_d+d1, i2*new_d+d2) += \
             GLG2_tmp(i1*new_d+d1, i2*GammaB.dim_d2+d)*YTiU(d, d2);
          }
        }
      } 
    }
  }
  //Check if new GLG2 matrix is really in canonical form:
  { 
    Eigen::VectorXcd id=Eigen::VectorXcd::Zero(new_d*new_d);
    for(int d=0; d<new_d; d++){  
      id(d*new_d+d)=1.;
    }
    id.normalize();
    Eigen::VectorXcd v2=PT_iTEBD_R_times(LambdaB, GLG2, id);
    std::complex<double> ritz=id.dot(v2);
    std::cout<<"Checking canonical form: ritz value="<<ritz<<" |v2/ritz-id|="<<(v2/ritz-id).norm()<<std::endl;
    v2=PT_iTEBD_LT_times(LambdaB, GLG2, id);
    ritz=id.dot(v2);
    std::cout<<"Checking canonical form: ritz value="<<ritz<<" |v2/ritz-id|="<<(v2/ritz-id).norm()<<std::endl;
  }
 }
  

  //Set sigma matrix
  int new_d=LambdaB.size(); 
  Eigen::MatrixXcd Sigma(dim_i*new_d, dim_i*new_d);
  for(int i1=0; i1<dim_i; i1++){
    for(int i2=0; i2<dim_i; i2++){
      for(int d1=0; d1<new_d; d1++){
        for(int d2=0; d2<new_d; d2++){
           Sigma(i1*new_d+d1, i2*new_d+d2) = LambdaB(d1)*GLG2(i1*new_d+d1, i2*new_d+d2)*LambdaB(d2);
        }
      }
    }
  }



  //Split Sigma matrix   
  time_point SVD2_start=now();
  TruncatedSVD_Return ret2=trunc.compress(Sigma);
  std::cout<<"second SVD: runtime="<<time_diff(now()-SVD2_start)<<"ms"<<std::endl;
  int new_d2=ret2.sigma.size();
  LambdaA=ret2.sigma;

  
  //Set new Gammas:
  GammaA.resize(dim_i, new_d, new_d2);
  for(int i=0; i<dim_i; i++){
    for(int d1=0; d1<new_d; d1++){
      for(int d2=0; d2<new_d2; d2++){
        GammaA(i,d1,d2)=ret2.U(i*new_d+d1, d2)*LambdaB(0)/LambdaB(d1);
      }
    }
  }
  GammaB.resize(dim_i, new_d2, new_d);
  for(int i=0; i<dim_i; i++){
    for(int d1=0; d1<new_d2; d1++){
      for(int d2=0; d2<new_d; d2++){
        GammaB(i,d1,d2)=ret2.Vdagger(d1, i*new_d+d2)*LambdaB(0)/LambdaB(d2);
      }
    }
  }
  LambdaB/=LambdaB(0)*LambdaB(0);

  //Balance weight between LambdaA and LambdaB:
  std::cout<<"GammaA.norm()="<<GammaA.norm()<<" GammaB.norm()="<<GammaB.norm()<<" LambdaA.norm()="<<LambdaA.norm()<<" LambdaB.norm()="<<LambdaB.norm()<<std::endl;
  std::complex<double> LambdaMean=sqrt(LambdaA(0)*LambdaB(0));
  LambdaA*=LambdaMean/LambdaA(0);
  LambdaB*=LambdaMean/LambdaB(0);

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


std::shared_ptr<ProcessTensorForward> PT_infinite(Parameters &param, DiagBB &diagBB){
  
  double dict_zero=param.get_as_double("dict_zero",0);
  TruncationLayout trunc_layout(param);


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
 
  infinite_normalize_specs specs;
  specs.iter=param.get_as_int("infinite_normalize_iter", 0); // 0: brute-force normalize, >0 power iteration, <0 Arnoldi
  specs.dont_normalize=param.get_as_bool("infinite_normalize_dont",false);
  specs.use_iter=param.get_as_bool("infinite_normalize_use_iter", specs.iter!=0);
  specs.use_Arnoldi=param.get_as_bool("infinite_normalize_Arnoldi", specs.use_iter);
  specs.eps=param.get_as_double("infinite_normalize_eps", 1e-15);

  Eigen::VectorXcd LambdaA(1); LambdaA(0)=1;
  Eigen::VectorXcd LambdaB(1); LambdaB(0)=1;
  MPS_Matrix GammaA(NL+1, 1,1); GammaA.fill(1.);
  MPS_Matrix GammaB(NL+1, 1,1); GammaB.fill(1.);

  for(int n=n_mem-1; n>0; n--){ 
    TruncatedSVD trunc=trunc_layout.get_base_line(n_mem-1-n, n_mem-1);
    trunc.keep=-1.;
    std::cout<<"------------------------------------------------------"<<std::endl;
    std::cout<<"n="<<n<<"/"<<n_mem<<" thr="<<trunc.threshold<<" dim_A="<<LambdaA.size()<<" dim_B="<<LambdaB.size()<<std::endl;
    std::cout<<"------------------------------------------------------"<<std::endl;
    Eigen::MatrixXcd expS=Eigen::MatrixXcd::Ones(NL+1,NL+1);
    expS.block(0,0,NL,NL)=diagBB.calculate_expS(n,tgrid.dt);
    if(n%2==0){
      if(specs.dont_normalize){
        PT_TEBD_step(LambdaA, LambdaB, GammaA, GammaB, expS, trunc);
      }else{
        PT_iTEBD_step(LambdaA, LambdaB, GammaA, GammaB, expS, trunc, specs);
      }
    }else{
      if(specs.dont_normalize){
        PT_TEBD_step(LambdaB, LambdaA, GammaB, GammaA, expS, trunc);
      }else{
        PT_iTEBD_step(LambdaB, LambdaA, GammaB, GammaA, expS, trunc, specs);
      }
    }
  }

  //First time step;
  std::cout<<"------------------------------------------------------"<<std::endl;
  std::cout<<"n="<<0<<"/"<<n_mem<<" dim_A="<<LambdaA.size()<<" dim_B="<<LambdaB.size()<<std::endl;
  std::cout<<"------------------------------------------------------"<<std::endl;
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


  //write:
  std::string write_PT=param.get_as_string("write_PT");
  int buffer_blocksize=param.get_as_int("buffer_blocksize",-1);

  if(param.get_as_bool("write_as_PTB",false)){
    std::shared_ptr<ProcessTensorForward> result(new ProcessTensorBuffer());
    ProcessTensorBuffer *PTB = dynamic_cast<ProcessTensorBuffer*>(result.get());
    if(write_PT!=""){
      PTB->set_new_file(write_PT, buffer_blocksize);
    }
    PTB->resize(tgrid.n_tot);
    for(int n=0; n<tgrid.n_tot; n++){
      ProcessTensorElement &e = PTB->get(n);
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
    PTB->get(0).M.inner_multiply_left(vl.transpose());
    PTB->get(tgrid.n_tot-1).M.inner_multiply_right(vr);
    PTB->expand_DiagBB(diagBB, dict_zero);
    PTB->calculate_closures();
    
    return result;
  }else{ //write as ProcessTensorRepeat
    std::shared_ptr<ProcessTensorForward> result(new ProcessTensorRepeat());
    ProcessTensorRepeat *PTR = dynamic_cast<ProcessTensorRepeat*>(result.get());
    PTR->set_specs(write_PT, buffer_blocksize);

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
  
    PTR->initial.resize(1);
    PTR->initial.get(0)=PTB.get(0);
    PTR->repeated.resize(1);
    PTR->repeated.get(0)=PTB.get(1);
   
    return result;
  }

}

}//namespace
