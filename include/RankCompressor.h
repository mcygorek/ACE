#ifndef RANK_COMPRESSOR_DEFINED_H
#define RANK_COMPRESSOR_DEFINED_H

#include <ctime>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/SVD>
#include "RRQR.h"
#include "MPS_Matrix.h"
#include "otimes.h"
#include "Parameters.h"

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
  double sum_threshold;
  int maxk;
  int reorthogonalize;
  int count, DUMP_SVD;
  bool use_BDCSVD;

  //increase threshold logarithmically from 'threshold' to 'threshold_to' within 'Nrange' SVDs
  double threshold_to;
  int Nrange;
 

  // check if parameters are set up so that any compression can happen
  bool has_effect()const{
    if(threshold>0. || sum_threshold>0. || maxk>0)return true;
    else return false;
  }

  //choose how many singular values should be kept:  
  int get_new_dim(const Eigen::VectorXd &sv){

    int newdim=sv.size();
    if(sv.size()<1){
      std::cerr<<"Compress (SVD): svd.singularValues().size()<1!"<<std::endl;
      exit(1);
    }

    // truncate based on given maximal inner dimension
    if(maxk>0 && maxk<newdim)newdim=maxk;

    // truncate based on (relative) magnitude of singular values:
    double this_threshold=threshold;
    if(Nrange>0){
      if(count<=0){
        this_threshold=threshold;
      }else if(count<=Nrange){
        this_threshold=threshold_to;
      }else{
        double x=(double)count/(double)Nrange;
        this_threshold=exp(log(threshold)*(1.-x)+log(threshold_to)*x);
std::cout<<"this threshold: "<<this_threshold<<std::endl;
      }
    }
    
    double cutoff=sv(0)*this_threshold;
    for(int i=1; i<newdim; i++){
      if(sv(i)<cutoff){newdim=i; break;}
    }

    // truncate based on (relative) sum of neglected SVs:
    if(sum_threshold>0){ 
      if(sum_threshold>1){
        std::cerr<<"Compress (SVD): sum_threshold>1!"<<std::endl;
        exit(1);
      }
      double sum=0.;
      for(int i=sv.size()-1; i>=0; i--){
        sum+=sv(i);
      }
      double cutoff=sum*sum_threshold;
      double newsum=0.;
      for(int i=newdim-1; i>=0; i--){
        newsum+=sv(i);
        if(newsum>=cutoff){
          newdim=i+1;
          break;
        }
      }
    }


    return newdim;
  }
   
  template <typename T>
  void compress_template(Eigen::MatrixXcd &A, Eigen::MatrixXcd &L, Eigen::MatrixXcd &R, bool low_to_high){

//    Eigen::JacobiSVD<Eigen::MatrixXcd>
    T svd( A, Eigen::ComputeThinU | Eigen::ComputeThinV);

    int newdim=get_new_dim(svd.singularValues());


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


    if(DUMP_SVD>=0 && count==DUMP_SVD){
//Note: Total number of compressions: (2*nmax-3)*N -> middle of chain for last step, backward direction: count ~ (2*nmax-3)*N-nmax/2-1  or (2*nmax-3)*(N-1/4)
      std::cout<<"Dumping SVD at "<<count<<"-th compression to file 'DUMP_SVD.dat'!"<<std::endl;
      std::ofstream ofs("DUMP_SVD.dat");
      for(int i=0; i<svd.singularValues().size(); i++){
        ofs<<svd.singularValues()(i)<<std::endl;
      }
    }

#ifdef PRINT_SVD_MATRIX_DIFF
std::cout<<"A-U*s*V^+:"<<std::endl;
std::cout<<A - svd.matrixU().block(0,0,A.rows(), newdim)*
               svd.singularValues().head(newdim).asDiagonal()*
               svd.matrixV().block(0,0,A.cols(),newdim).adjoint()<<std::endl;
#endif




#ifdef PRINT_SVD_ORTHO
{
std::ofstream ofs("PRINT_SVD_ORTHO", std::ios_base::app);
std::time_t print_svd_ortho_time=std::time(NULL);
double diffU=max_diff_from_ortho(svd.matrixU());
double diffV=max_diff_from_ortho(svd.matrixV());
ofs<<diffU<<" "<<diffV<<" "<<std::asctime(std::localtime(&print_svd_ortho_time))<<std::endl;
/*
if(diffU>1e-12){
std::cout<<"SVD produces non-orthogonal U: "<<diffU<<std::endl;
print_diff_from_ortho(svd.matrixU());
}
if(diffV>1e-12){
std::cout<<"SVD produces non-orthogonal V: "<<diffV<<std::endl;
print_diff_from_ortho(svd.matrixV());
}
*/
if(diffU>1e-12||diffV>1e-12){
std::cout<<"newdim="<<newdim<<std::endl;
std::cout<<"A-U*s*V^+:"<<std::endl;
Eigen::MatrixXcd diffmat = A - svd.matrixU().block(0,0,A.rows(), newdim)*
                           svd.singularValues().head(newdim).asDiagonal()*
                           svd.matrixV().block(0,0,A.cols(),newdim).adjoint();
{std::ofstream ofs2("PROBLEMATIC_MATRIX.dat"); ofs2<<A;}
{std::ofstream ofs2("PROBLEMATIC_SVs.dat"); ofs2<<svd.singularValues();}
double maxelem=0;
int maxi=0, maxj=0;
for(int i=0; i<diffmat.rows(); i++){
  for(int j=0; j<diffmat.cols(); j++){
    if(abs(diffmat(i,j))>maxelem){
      maxelem=abs(diffmat(i,j));
      maxi=i; maxj=j;
    }
  }
}
std::cout<<maxi<<" "<<maxj<<" "<<diffmat(maxi,maxj)<<std::endl;
exit(1);
}
}
#endif

    if(reorthogonalize>0){
      if(low_to_high){
        Eigen::MatrixXcd U=svd.matrixU().block(0,0,A.rows(),newdim);
        for(int run=0; run<reorthogonalize; run++){
          for(int i=0; i<U.cols(); i++){
            for(int j=0; j<i; j++){
              U.col(i)-=(U.col(j).dot(U.col(i)))*U.col(j);
            }
            U.col(i).normalize();
          }
        }
        L = U;
        R = U.adjoint() * A;
      }else{
        Eigen::MatrixXcd V=svd.matrixV().block(0,0,A.cols(),newdim);
        for(int run=0; run<reorthogonalize; run++){
          for(int i=0; i<V.cols(); i++){
            for(int j=0; j<i; j++){
              V.col(i)-=(V.col(j).dot(V.col(i)))*V.col(j);
            }
            V.col(i).normalize();
          }
        }
        R = V.adjoint();
        L = A * V;
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

    count++;
  }

  virtual void compress(Eigen::MatrixXcd &A, Eigen::MatrixXcd &L, Eigen::MatrixXcd &R, bool low_to_high){

    if(use_BDCSVD){
      compress_template<Eigen::BDCSVD<Eigen::MatrixXcd> >(A, L, R, low_to_high);
    }else{
      compress_template<Eigen::JacobiSVD<Eigen::MatrixXcd> >(A, L, R, low_to_high);
    }
  }
 

  void setup(Parameters &param){
    threshold=param.get_as_double("threshold", 0);
    sum_threshold=param.get_as_double("sum_threshold", 0);
    maxk=param.get_as_size_t("compress_maxk", 0);
    reorthogonalize=param.get_as_int("reorthogonalize",0);
    DUMP_SVD=param.get_as_int("compress_dump",-1);
    count=param.get_as_size_t("compress_dump_initial_count",0);
    use_BDCSVD=param.get_as_bool("use_BDCSVD",false);

    if(param.is_specified("threshold_range")){
      threshold=param.get_as_double("threshold_range",0.,0,0);
      threshold_to=param.get_as_double("threshold_range",0.,0,1);
      Nrange=param.get_as_int("threshold_range",0.,0,2);
    }else{
      Nrange=0;
    }
  }
  RankCompressor_SVD(Parameters &param){
    setup(param);
  }
  RankCompressor_SVD(double thresh=0., double sum_thr_=0., int maxk_=0, int reortho=0, int dump=-1)
   : threshold(thresh),sum_threshold(sum_thr_), maxk(maxk_),reorthogonalize(reortho), count(0), DUMP_SVD(dump){
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
