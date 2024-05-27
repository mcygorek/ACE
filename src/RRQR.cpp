#include "RRQR.hpp"
#include <Eigen/Core>
#include <Eigen/SVD>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>

namespace ACE{

void HouseholderMultiply(const Eigen::VectorXcd &v,Eigen::MatrixXcd &A, int start_col){
//std::cout<<"HH: "<<v.size()<<" "<<A.rows()<<" "<<A.cols()<<std::endl;
//  Eigen::VectorXcd b=v.adjoint()*A;
  Eigen::VectorXcd b;
  if(v.size()==A.rows()){
    b=v.adjoint()*A;
  }else{
    b=v.head(A.rows()).adjoint()*A;
  }
  for(int i=0; i<A.rows(); i++){
    for(int j=start_col; j<A.cols(); j++){
      A(i,j)-=2.*v(i)*b(j);
    }
  }
}

std::vector<int> get_reverse_pivot(const std::vector<int> &pivot, bool print_pivots){
  //reverse pivot lookup
  std::vector<int> revpiv(pivot.size());
  for(size_t i=0; i<pivot.size(); i++){
    revpiv[pivot[i]]=i;
  }
  if(print_pivots){
    std::cout<<"pivots:"<<std::endl;
    for(size_t i=0; i<pivot.size(); i++){
      std::cout<<pivot[i]<<" ";
    }std::cout<<std::endl;
    std::cout<<"reverse pivots:"<<std::endl;
    for(size_t i=0; i<revpiv.size(); i++){
      std::cout<<revpiv[i]<<" ";
    }std::cout<<std::endl;
    std::cout<<std::endl;
  }
  return revpiv;
}

/** performs a rank-reducing factorization of a tall matrix A into L and R
    A is modified in the process */
void QRP(Eigen::MatrixXcd &A, Eigen::MatrixXcd &Q, Eigen::MatrixXcd &R, std::vector<int> &pivot, double threshold, int debug_mode){

  int maxstep=A.cols();
  if(A.rows()<A.cols()){
    maxstep=A.rows();
  }
  if(A.cols()<1||A.rows()<1){
    Q.resize(A.rows(), A.cols());
    R.resize(A.cols(), A.cols());
    return;
  }

/*
  if(A.cols()==1){
//    double n=A.col(0).norm();
//    if(n<1e-16)n=1.;
    double n=1.;
    Q=A/n;
    R.resize(1,1); R(0,0)=n;
    return; 
  }
*/

//  bool silent=!(debug_mode&1);
  bool nopivot=(debug_mode &16);
  int threshold_reached=-1;
//  double threshold_remainder=0.;
  
//-----------------------
  
  //storage for column pivots:
  pivot.resize(A.cols());
  for(size_t i=0; i<pivot.size(); i++)pivot[i]=i;

  //storage for Householder vectors
  std::vector<Eigen::VectorXcd> V;

  //storage of column norm squares (for deciding on pivoting)

  std::vector<double> sn(A.cols(),0.);
  for(int i=0; i<A.cols(); i++){
    for(int j=0; j<A.rows(); j++){
      sn[i]+=(A(j,i).real()*A(j,i).real()+A(j,i).imag()*A(j,i).imag());
    }
  }

  double cutoff=threshold;
  int rank=maxstep;
  for(int step=0; step<maxstep; step++){
//    if(!silent)std::cout<<"step: "<<step<<std::endl;

    //pivot
    if(!nopivot){
      //determine pivot 
      int m=step;
      for(int j=step+1; j<A.cols(); j++){
        if(sn[j]>sn[m])m=j;
      }

      //apply pivot
      if(m!=step){
        int tmp=pivot[step];
        pivot[step]=pivot[m];
        pivot[m]=tmp;

        Eigen::VectorXcd swp=A.col(step);
        A.col(step)=A.col(m);
        A.col(m)=swp;
    
        double tmp_=sn[step];
        sn[step]=sn[m];
        sn[m]=tmp_;
      }
    }

    //determine Householder vector
    Eigen::VectorXcd v=A.col(step);
    for(int i=0; i<step; i++)v(i)=0;
//    std::complex<double> lambda=-v.norm()*(v(step)/abs(v(step)));
    std::complex<double> lambda;
    if(v(step).real()*v(step).real()+v(step).imag()*v(step).imag()>1e-20){
      lambda=-v.norm()*(v(step)/abs(v(step)));
    }else{
      lambda=-v.norm();
    }

    v(step)-=lambda;
    v.normalize();
    V.push_back(v);
    HouseholderMultiply(v,A,step);

    //adjust column moduli
    for(int j=step+1; j<A.cols(); j++){
      sn[j]-=(A(step,j).real()*A(step,j).real()+A(step,j).imag()*A(step,j).imag());
/*
      double d=0; 
      for(int l=step+1; l<A.rows(); l++){
        d+=(A(l,j).real()*A(l,j).real()+A(l,j).imag()*A(l,j).imag());
      }
      if(fabs(sn[j]-d)>1e-8)std::cout<<"sn["<<j<<"] "<<sn[j]<<" "<<d<<" "<<sn[j]-d<<std::endl;*/
    }


    //stop if threshold has been reached
    if(step==0){
      cutoff=abs(A(0,0))*threshold;
    }else if(threshold_reached<0 && abs(A(step,step))<cutoff){ //threshold){
      threshold_reached=step;
      rank=threshold_reached+1;
//      if(!silent)std::cout<<"Break at step: "<<step<<std::endl;
      break;
    }
  }
//  if(!silent)std::cout<<"rank: "<<rank<<std::endl;


  //Build final Q and R matrices
  Q=Eigen::MatrixXcd::Identity(A.rows(), rank);
  for(int i=(int)V.size()-1; i>=0; i--){
    HouseholderMultiply(V[i], Q);
  }


  R=Eigen::MatrixXcd(rank, A.cols());
  for(int i=0; i<rank; i++){
    for(int j=0; j<A.cols(); j++){
      R(i,j)=A(i,j);
//      R(i,j)=A(i,revpiv[j]);
    }
  }
}

void QLV(Eigen::MatrixXcd &A, Eigen::MatrixXcd &Q, Eigen::MatrixXcd &L, Eigen::MatrixXcd &V,double threshold, int debug_mode){


  std::vector<int> pivot1;
  Eigen::MatrixXcd Q1, Rtmp1;
  QRP(A, Q1, Rtmp1, pivot1, threshold*0.01, debug_mode);
  std::vector<int> revpiv1=get_reverse_pivot(pivot1);
  
//std::cout<<"Intermediate dim: "<<Rtmp1.rows()<<std::endl;
  Eigen::MatrixXcd Q2, Rtmp2;
  std::vector<int> pivot2;
  Eigen::MatrixXcd Rtmp_T=Rtmp1.transpose();
  QRP(Rtmp_T, Q2, Rtmp2, pivot2, threshold, debug_mode); //| 16);

  //R^T= Q_2 * R_2 * P_2^T -> R = P_2 * R_2^T * Q_2^T
  //A = Q R P^T -> Q*P_2 * R_2^T * Q_2^T*P^T

  Q=Eigen::MatrixXcd(Q1.rows(), Q1.cols());
  for(int i=0; i<Q1.rows(); i++){
    for(int j=0; j<Q1.cols(); j++){
      Q(i,j)=Q1(i,pivot2[j]);
    }
  }
  V=Eigen::MatrixXcd(Q2.cols(), Q2.rows());
  for(int i=0; i<Q2.rows(); i++){
    for(int j=0; j<Q2.cols(); j++){
      V(j,i)=Q2(revpiv1[i],j);
    }
  }
  L=Rtmp2.transpose();
}

void QLV_with_debug(Eigen::MatrixXcd &A, Eigen::MatrixXcd &Q, Eigen::MatrixXcd &L, Eigen::MatrixXcd &V,double threshold, int debug_mode, int times){

  bool print_matrices=(debug_mode&8);
  bool test_result=(debug_mode&2);
//  bool print_pivots=(debug_mode&4);
  bool compare_svd=(debug_mode &32);
  if(print_matrices){
    std::cout<<"Original matrix:"<<std::endl;
    std::cout<<A<<std::endl<<std::endl;
  }
  Eigen::MatrixXcd A_backup;
  if(test_result){
    A_backup=A;
  } 

//  Eigen::ColPivHouseholderQR<Eigen::MatrixXcd> qr(Rtmp_T);
//  Rtmp2=qr.matrixR().triangularView<Eigen::Upper>();
//  Q2=qr.householderQ();


  QLV(A, Q, L, V, threshold, debug_mode);

  for(int i=1; i<times; i++){
    Eigen::MatrixXcd Qtmp, Ltmp, Vtmp;
    QLV(L, Qtmp, Ltmp, Vtmp, threshold, debug_mode);
    Q=Q*Qtmp; V=Vtmp*V; L=Ltmp;
  }
 


//-----------------------

  if(print_matrices){
    std::cout<<"Q^dagger*Q:"<<std::endl<<Q.adjoint()*Q<<std::endl<<std::endl;
    std::cout<<"Q*L*V:"<<std::endl<<Q*L*V<<std::endl<<std::endl;
  }
#ifdef PRINT_SVD
    if(L.rows()<L.cols())std::cout<<"RRQR("<<L.rows()<<"):";
    else std::cout<<"RRQR("<<L.cols()<<"):";
    for(int i=0; i<L.rows() && i<L.cols(); i++){
      std::cout<<" "<<std::abs(L(i,i));
    }std::cout<<std::endl;
#endif
  if(compare_svd){
    std::ofstream ofs_diags("DUMP_R_DIAGS.dat");
    for(int i=0; i<L.rows() && i<L.cols(); i++){
      ofs_diags<<std::abs(L(i,i))<<std::endl;
    }
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(A_backup);
    std::ofstream ofs_svd("DUMP_SVD.dat");
    for(int i=0; i<svd.singularValues().size(); i++){
      ofs_svd<<svd.singularValues()(i)<<std::endl;
    }
  }

  if(test_result){
    Eigen::MatrixXcd B=Q*L*V;
    double maxdiff=0.;
    for(int i=0; i<A_backup.rows(); i++){
      for(int j=0; j<A_backup.cols(); j++){
        double a=abs(A_backup(i,j)-B(i,j));
        if(a>maxdiff)maxdiff=a;
      }
    }
    std::cout<<"Maximal difference between A and QR: "<<maxdiff<<std::endl<<std::endl;
    if(maxdiff>threshold*100){
      std::ofstream ofs("DUMP_A.dat");
      for(int i=0; i<A_backup.rows(); i++){
        for(int j=0; j<A_backup.cols(); j++){
          ofs<<A_backup(i,j)<<" ";
        }
        ofs<<std::endl;
      }
      exit(1);
    }
  }
}


void RRQR(Eigen::MatrixXcd &A, Eigen::MatrixXcd &L, Eigen::MatrixXcd &R, double threshold, int debug_mode, int times){

  Eigen::MatrixXcd Q, Ltmp, V;
  QLV_with_debug(A, Q, Ltmp, V, threshold, debug_mode, times);
  L=Q;
  R=Ltmp*V;

}

}//namespace
