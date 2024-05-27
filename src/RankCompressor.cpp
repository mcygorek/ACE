#include "RankCompressor.hpp"
#include <ctime>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/SVD>
#include "RRQR.hpp"
#include "MPS.hpp"
#include "otimes.hpp"
#include "Parameters.hpp"
#include "DummyException.hpp"

namespace ACE{

  //Note: left_to_right=true means: start from a[0]
template <typename T>
  void RankCompressor_ScalarType<T>::compress(MPS_Matrix_ScalarType<T> &a, 
               Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &L, 
               Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &R, 
               bool low_to_high, double s){

    if(fabs(s)<1e-30)s=sqrt(sqrt(a.dim_i));

    if(low_to_high){
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A(a.dim_d1*a.dim_i, a.dim_d2);
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
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A(a.dim_d1, a.dim_d2*a.dim_i);
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

template <typename T>
  void RankCompressor_ScalarType<T>::compress_keep_largest(
                     MPS_Matrix_ScalarType<T> &a, 
                     Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &L, 
                     Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &R,   
                     bool low_to_high){

    if(low_to_high){
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A(a.dim_d1*a.dim_i, a.dim_d2);
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
        throw DummyException();
      }
      L*=R(0,0);
      R/=R(0,0);
    }else{
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A(a.dim_d1, a.dim_d2*a.dim_i);
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
        throw DummyException();
      }
      L/=L(0,0);
      R*=L(0,0);
    }
  }
 
template <typename T>
  void RankCompressor_ScalarType<T>::sweep_block_low_to_high(
                int n, MPS_ScalarType<T> &mps, double keep_weight, 
                Sweep_Trafo_Processor_ScalarType<T> *proc){

#ifdef PRINT_SVD
std::cout<<"sweep low->high: "<<n<<std::endl;
#endif

    if(n<0)return;
    if(n>=(int)mps.a.size()-1)return;

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> L, R;
    compress(mps.a[n], L, R, true, keep_weight);
    int newdim=L.cols();

    //modify a[n]
    {
    MPS_Matrix_ScalarType<T> ar(mps.a[n].dim_i, mps.a[n].dim_d1, newdim);
    for(int d1=0; d1<mps.a[n].dim_d1; d1++){
      for(int i=0; i<mps.a[n].dim_i; i++){
        for(int dn=0; dn<newdim; dn++){
          ar(i,d1,dn) = L(d1*mps.a[n].dim_i+i, dn);
        }
      }
    }
    mps.a[n].swap(ar);   
    }
   
    //modify a[n+1]
    {
    MPS_Matrix_ScalarType<T> ar(mps.a[n+1].dim_i, newdim, mps.a[n+1].dim_d2);
    ar.set_zero();
    for(int d1=0; d1<mps.a[n+1].dim_d1; d1++){
      for(int d2=0; d2<mps.a[n+1].dim_d2; d2++){
        for(int i=0; i<mps.a[n+1].dim_i; i++){
          for(int dn=0; dn<newdim; dn++){
            ar(i,dn,d2) += R(dn, d1) * mps.a[n+1](i, d1, d2);
          }
        }
      }
    }
    mps.a[n+1].swap(ar);   
    }
   
    //apply to environment operators
    if(proc!=NULL){ 
      proc->process_low_to_high(n, R);
    }
  }

template <typename T>
  void RankCompressor_ScalarType<T>::sweep_block_high_to_low(
            int n, MPS_ScalarType<T> &mps, double keep_weight, 
            Sweep_Trafo_Processor_ScalarType<T> *proc){

    if(n<1)return;
#ifdef PRINT_SVD
std::cout<<"sweep high->low: "<<n<<std::endl;
#endif

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> L, R;

    compress(mps.a[n], L, R, false, keep_weight);

    int newdim=L.cols();

    //modify a[n]
    {
    MPS_Matrix_ScalarType<T> ar(mps.a[n].dim_i, newdim, mps.a[n].dim_d2);
    for(int d2=0; d2<mps.a[n].dim_d2; d2++){
      for(int i=0; i<mps.a[n].dim_i; i++){
        for(int dn=0; dn<newdim; dn++){
          ar(i,dn,d2) = R(dn, i*mps.a[n].dim_d2+d2);
        }
      }
    }
    mps.a[n].swap(ar);   
    }
   
    //modify a[n-1]
    {
    MPS_Matrix_ScalarType<T> ar(mps.a[n-1].dim_i, mps.a[n-1].dim_d1, newdim);
    ar.set_zero();
    for(int d1=0; d1<mps.a[n-1].dim_d1; d1++){
      for(int d2=0; d2<mps.a[n-1].dim_d2; d2++){
        for(int i=0; i<mps.a[n-1].dim_i; i++){
          for(int dn=0; dn<newdim; dn++){
            ar(i,d1,dn) += mps.a[n-1](i, d1, d2) * L(d2, dn);
          }
        }
      }
    }
    mps.a[n-1].swap(ar);   
    }

    //apply to environment operators
    if(proc!=NULL){ 
      proc->process_high_to_low(n, L);
    }
  }

template <typename T>
  void RankCompressor_ScalarType<T>::low_end_multiply_and_compress(
      MPS_ScalarType<T> &mps, const MPS_ScalarType<T> & other, int sweep_start){

    if(mps.a.size()<other.a.size()){
      std::cerr<<"MPS::multiply_front: a.size()<other.a.size()!"<<std::endl;
      throw DummyException();
    }
    for(size_t i=0; i<other.a.size(); i++){
std::cout<<"forward sweep "<<i<<"/"<<other.a.size()<<std::endl;

      mps.a[i].multiply(other.a[i]);

      if((int)i>sweep_start){
        bool skip=false; 
        if(mps.a[i-1].dim_d1==1 && mps.a[i-1].dim_d2==1)skip=true;
        if(other.a[i-1].dim_d1==1 && other.a[i-1].dim_d2==1)skip=true;
        if(!skip){
          sweep_block_low_to_high(i-1, mps);
        }
      }
    }
    for(int i=other.a.size()-1; i>1; i--){
std::cout<<"backward sweep "<<i<<"/"<<other.a.size()<<std::endl;
      sweep_block_high_to_low(i, mps);
    }
  }

template <typename T>
  void RankCompressor_ScalarType<T>::high_end_multiply_and_compress(
            MPS_ScalarType<T> & mps, const MPS_ScalarType<T> & other, 
            double keep_weight, Compress_Trafo_At *cta){

    if(mps.a.size()<other.a.size()){
      std::cerr<<"MPS::high_end_multiply_and_compress: a.size()<other.a.size()!"<<std::endl;
      throw DummyException();
    }

    int start=((int)mps.a.size())-((int)other.a.size());
    
    for(int i=((int)other.a.size())-1; i>=0; i--){
      mps.a[i+start].multiply(other.a[i]);

      if(i<((int)other.a.size())-1){ 
        if(cta!=NULL && i+start+1==cta->n+1){
          std::cerr<<"MPS::high_end_multiply_and_compress: using cta not implemented yet!"<<std::endl;
          throw DummyException();
//          sweep_block_high_to_low(i+start+1, mps, keep_weight, &cta->L);
        }else{
          sweep_block_high_to_low(i+start+1, mps, keep_weight);
        } 
      }
    }
    for(int i=start; i<(int)mps.a.size()-1; i++){
      if(cta!=NULL && i==cta->n){
        std::cerr<<"MPS::high_end_multiply_and_compress: using cta not implemented yet!"<<std::endl;
        throw DummyException();
//        sweep_block_low_to_high(i, mps, keep_weight, &cta->R);
      }else{
        sweep_block_low_to_high(i, mps, keep_weight);
      }
    }
  }

template <typename T>
  void RankCompressor_ScalarType<T>::multiply_and_compress(
               MPS_ScalarType<T> & mps, const MPS_ScalarType<T> & other){

    if(mps.a.size()!=other.a.size()){
      std::cerr<<"MPS::multiply_front: a.size()<other.a.size()!"<<std::endl;
      throw DummyException();
    }
//    low_end_multiply_and_compress(other, mps);
    low_end_multiply_and_compress(mps, other);
  }

template <typename T>
  void RankCompressor_ScalarType<T>::sweep(MPS_ScalarType<T> & mps, double keep_weight){
     for(int n=0; n<(int)mps.a.size()-1; n++){
       sweep_block_low_to_high(n, mps, keep_weight);
     }
     for(int n=(int)mps.a.size()-1; n>1; n--){
       sweep_block_high_to_low(n, mps, keep_weight);
     }
  }


template class RankCompressor_ScalarType<std::complex<double>>;
template class RankCompressor_ScalarType<double>;

template <typename T>
void RankCompressor_None_ScalarType<T>::compress(
                        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A, 
                        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &L, 
                        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &R, 
                        bool low_to_high){
    if(low_to_high){
      L=A;
      R=Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(A.cols(), A.cols());
    }else{
      L=Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(A.rows(), A.rows());
      R=A;
    }
  }

template class RankCompressor_None_ScalarType<std::complex<double>>;
template class RankCompressor_None_ScalarType<double>;

}//namespace
