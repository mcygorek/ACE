#ifndef RANK_COMPRESSOR_DEFINED_H
#define RANK_COMPRESSOR_DEFINED_H

#include <cmath>
#include <Eigen/Core>
#include <Eigen/SVD>
#include "RRQR.h"
#include "MPS_Matrix.h"

class RankCompressor{
public:
  int debug;
  bool keep_weight;

  //Note: left_to_right=true means: start from a[0]
  virtual void compress(Eigen::MatrixXcd &A, Eigen::MatrixXcd &L, Eigen::MatrixXcd &R, bool low_to_high)=0;

  virtual void compress(MPS_Matrix &a, Eigen::MatrixXcd &L, Eigen::MatrixXcd &R, bool low_to_high, double s=0.){

    if(fabs(s)<1e-30)s=sqrt(sqrt(a.dim_i));

    if(low_to_high){
      Eigen::MatrixXcd A(a.dim_d1*a.dim_i, a.dim_d2);
      for(int d1=0; d1<a.dim_d1; d1++){
        for(int i=0; i<a.dim_i; i++){
          for(int d2=0; d2<a.dim_d2; d2++){
            A(d1*a.dim_i+i, d2) = a(i, d1, d2);
          }
        }
      }
      compress(A, L, R, low_to_high);
      if(keep_weight && R.rows()>0){
//        double s=std::abs(R(0,0));
//std::cout<<std::endl<<"s: "<<s<<std::endl<<std::endl;
        if(s>1.){
//          s=a.dim_i;
          L*=s;
          R/=s;
        }
      }
    }else{
      Eigen::MatrixXcd A(a.dim_d1, a.dim_d2*a.dim_i);
      for(int d1=0; d1<a.dim_d1; d1++){
        for(int i=0; i<a.dim_i; i++){
          for(int d2=0; d2<a.dim_d2; d2++){
            A(d1,i*a.dim_d2+d2) = a(i, d1, d2);
          }
        }
      }
      compress(A, L, R, low_to_high);

      if(keep_weight && L.cols()>0){
//        double s=std::abs(L(0,0));
        if(s>1){
//          s=a.dim_i;
          L/=s;
          R*=s;
        }
      }
    }
  }
  virtual void compress_keep_largest(MPS_Matrix &a, Eigen::MatrixXcd &L, Eigen::MatrixXcd &R, bool low_to_high){

    if(low_to_high){
      Eigen::MatrixXcd A(a.dim_d1*a.dim_i, a.dim_d2);
      for(int d1=0; d1<a.dim_d1; d1++){
        for(int i=0; i<a.dim_i; i++){
          for(int d2=0; d2<a.dim_d2; d2++){
            A(d1*a.dim_i+i, d2) = a(i, d1, d2);
          }
        }
      }
      compress(A, L, R, low_to_high);
      if(R.rows()<1 || R.cols()<1){
        std::cerr<<"compress_keep_largest: R.rows()<1 || R.cols()<1!"<<std::endl;
        exit(1);
      }
      L*=R(0,0);
      R/=R(0,0);
    }else{
      Eigen::MatrixXcd A(a.dim_d1, a.dim_d2*a.dim_i);
      for(int d1=0; d1<a.dim_d1; d1++){
        for(int i=0; i<a.dim_i; i++){
          for(int d2=0; d2<a.dim_d2; d2++){
            A(d1,i*a.dim_d2+d2) = a(i, d1, d2);
          }
        }
      }
      compress(A, L, R, low_to_high);
      if(L.rows()<1 || L.cols()<1){
        std::cerr<<"compress_keep_largest: L.rows()<1 || L.cols()<1!"<<std::endl;
        exit(1);
      }
      L/=L(0,0);
      R*=L(0,0);
    }
  }
 
  RankCompressor(){
    debug=0;
    keep_weight=true;
  }
  virtual ~RankCompressor(){
  }
};


class RankCompressor_None: public RankCompressor{
public:
  virtual void compress(Eigen::MatrixXcd &A, Eigen::MatrixXcd &L, Eigen::MatrixXcd &R, bool low_to_high){
    if(low_to_high){
      L=A;
      R=Eigen::MatrixXcd::Identity(A.cols(), A.cols());
    }else{
      L=Eigen::MatrixXcd::Identity(A.rows(), A.rows());
      R=A;
    }
  }
};

class RankCompressor_SVD: public RankCompressor{
public:
  double threshold;
  int maxk;
  bool reorthogonalize;

    
  virtual void compress(Eigen::MatrixXcd &A, Eigen::MatrixXcd &L, Eigen::MatrixXcd &R, bool low_to_high){


    Eigen::JacobiSVD<Eigen::MatrixXcd> svd( A, Eigen::ComputeThinU | Eigen::ComputeThinV);

    int newdim=svd.singularValues().size();
    if(svd.singularValues().size()<1){
      std::cerr<<"Compress (SVD): svd.singularValues().size()<1!"<<std::endl;
      exit(1);
    }

    if(maxk>0 && maxk<newdim)newdim=maxk;

    double cutoff=svd.singularValues()(0)*threshold;
    for(int i=1; i<newdim; i++){
      if(svd.singularValues()(i)<cutoff){newdim=i; break;}
    }

#ifdef PRINT_SVD_DIMS
    std::cout<<"SVD dims: "<<A.rows()<<", "<<A.cols()<<" -> "<<newdim<<std::endl;
#endif
#ifdef PRINT_SVD
    std::cout<<"SVDs("<<newdim<<"):";
    for(int i=0; i<svd.singularValues().size(); i++){
      std::cout<<" "<<svd.singularValues()(i);
    }std::cout<<std::endl;
#elif defined(PRINT_SVD_FIRST)
   std::cout<<"SVDs("<<newdim<<"):";
   std::cout<<" first: "<<svd.singularValues()(0)<<std::endl;
#endif

#ifdef PRINT_SVD_MATRIX_DIFF
std::cout<<"A-U*s*V^+:"<<std::endl;
std::cout<<A - svd.matrixU().block(0,0,A.rows(), newdim)*
               svd.singularValues().head(newdim).asDiagonal()*
               svd.matrixV().block(0,0,A.cols(),newdim).adjoint()<<std::endl;
#endif


    if(reorthogonalize){
      Eigen::MatrixXcd V=svd.matrixV().block(0,0,A.cols(),newdim);

      if(low_to_high){
        Eigen::MatrixXcd U=svd.matrixU().block(0,0,A.rows(),newdim);
        for(int i=0; i<newdim; i++){
          for(int j=i-1; j>=0; j--){
            U.col(i)-=(U.col(j).dot(U.col(i)))*U.col(j);
          }
          U.col(i).normalize();
        }

        L = U;
        R = Eigen::MatrixXcd(newdim, A.cols());
        for(int i=0; i<newdim; i++){
          R.row(i)=U.col(i).adjoint()*A;
        } 

        if(true){
          for(int i=0; i<newdim; i++){
            double sn=sqrt(R.row(i).norm());
            R.row(i)/=sn;
            L.col(i)*=sn;
          }
//std::cout<<"."<<std::flush;
        }

      }else{
        Eigen::MatrixXcd V=svd.matrixV().block(0,0,A.cols(),newdim);
        for(int i=0; i<newdim; i++){
          for(int j=i-1; j>=0; j--){
            V.col(i)-=(V.col(j).dot(V.col(i)))*V.col(j);
          }
          V.col(i).normalize();
        }
 
        R = V.adjoint();
        L = Eigen::MatrixXcd(A.rows(), newdim);
        for(int i=0; i<newdim; i++){
          L.col(i)=A*V.col(i);
        }

        if(true){
          for(int i=0; i<newdim; i++){
            double sn=sqrt(L.col(i).norm());
            L.col(i)/=sn;
            R.row(i)*=sn;
          }
        }
      } 

    }else{
      if(low_to_high){
        L = svd.matrixU().block(0,0,A.rows(), newdim);
        R = svd.singularValues().head(newdim).asDiagonal() 
            * svd.matrixV().block(0,0,A.cols(),newdim).adjoint();
      }else{
        L = svd.matrixU().block(0,0,A.rows(), newdim)
            *svd.singularValues().head(newdim).asDiagonal() ;
        R = svd.matrixV().block(0,0,A.cols(),newdim).adjoint();
      }
    }

  }

  RankCompressor_SVD(double thresh, int maxk_=0, bool reortho=false)
   : threshold(thresh),maxk(maxk_),reorthogonalize(reortho){
  }
  virtual ~RankCompressor_SVD(){
  }
};



class RankCompressor_RRQR: public RankCompressor{
public: 
  double threshold;
  
  virtual void compress(Eigen::MatrixXcd &A, Eigen::MatrixXcd &L, Eigen::MatrixXcd &R, bool low_to_high){
    int times=1;
    Eigen::MatrixXcd Q, Ltmp, V;
    QLV_with_debug(A, Q, Ltmp, V, threshold, debug, times);

    if(low_to_high){
      L=Q;
      R=Ltmp*V;
    }else{
      L=Q*Ltmp;
      R=V;
    }
  }

  RankCompressor_RRQR(double thresh): threshold(thresh){
  }
  virtual ~RankCompressor_RRQR(){
  }
};






/*
class RankCompressor_SVD_Block: public RankCompressor_SVD{
public:
    
  virtual void compress(MPS_Matrix &a, Eigen::MatrixXcd &L, Eigen::MatrixXcd &R, bool left_to_right){

    L.resize(0,0);
    R.resize(0,0);
    std::vector<Eigen::MatrixXcd> Lvec(a.dim_i);
    std::vector<Eigen::MatrixXcd> Rvec(a.dim_i);
    if(left_to_right){
      int sum_Lrows=0, sum_Lcols=0, sum_Rcols;
      for(int i=0; i<a.dim_i; i++){
        Eigen::MatrixXcd A(a.dim_d1, a.dim_d2);
        for(int d1=0; d1<a.dim_d1; d1++){
          for(int d2=0; d2<a.dim_d2; d2++){
            A(d1, d2) = a(i, d1, d2);
          }
        }
        RankCompressor_SVD::compress(A, Lvec[i], Rvec[i], left_to_right);
      }
    }else{
      Eigen::MatrixXcd A(a.dim_d1, a.dim_d2*a.dim_i);
      for(int d1=0; d1<a.dim_d1; d1++){
        for(int i=0; i<a.dim_i; i++){
          for(int d2=0; d2<a.dim_d2; d2++){
            A(d1,i*a.dim_d2+d2) = a(i, d1, d2);
          }
        }
      }
      RankCompressor_SVD::compress(A, L, R, left_to_right);
    }
  }

  RankCompressor_SVD_Block(double thresh): RankCompressor_SVD(thresh){
  }
};
*/


#include "RankCompressor_PowerURV.h"

#endif
