#ifndef RANK_COMPRESSOR_SVD_DEFINED_H
#define RANK_COMPRESSOR_SVD_DEFINED_H

#include "RankCompressor.hpp"

namespace ACE{

template <typename T>
class RankCompressor_SVD_ScalarType: public RankCompressor_ScalarType<T>{
public:
  double threshold;
  double sum_threshold;
  int maxk;
  int reorthogonalize;
  int count, DUMP_SVD;
//  bool use_BDCSVD;

  //increase threshold logarithmically from 'threshold' to 'threshold_to' within 'Nrange' SVDs
  double threshold_to;
  int Nrange;
 

  // check if parameters are set up so that any compression can happen
  inline bool has_effect()const{
    if(threshold>0. || sum_threshold>0. || maxk>0)return true;
    else return false;
  }

  //choose how many singular values should be kept:  
  int get_new_dim(const Eigen::VectorXd &sv);

//  template <typename T2>
  void compress_template(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A, 
                         Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &L, 
                         Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &R, 
                         bool low_to_high);

  virtual void compress(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A, 
                        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &L, 
                        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &R, 
                        bool low_to_high);


  void setup(Parameters &param);

  inline RankCompressor_SVD_ScalarType(Parameters &param){
    setup(param);
  }
  inline RankCompressor_SVD_ScalarType(double thresh=0., double sum_thr_=0., int maxk_=0, int reortho=0, int dump=-1)
   : threshold(thresh),sum_threshold(sum_thr_), maxk(maxk_),reorthogonalize(reortho), count(0), DUMP_SVD(dump), Nrange(0){
  }
  inline virtual ~RankCompressor_SVD_ScalarType(){
  }
};

typedef RankCompressor_SVD_ScalarType<std::complex<double> > RankCompressor_SVD;
typedef RankCompressor_SVD_ScalarType<double> RankCompressor_SVD_real;


}//namespace
#endif
