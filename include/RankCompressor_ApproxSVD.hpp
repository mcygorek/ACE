#ifndef RANK_COMPRESSOR_APPROX_SVD_DEFINED_H
#define RANK_COMPRESSOR_APPROX_SVD_DEFINED_H

#include "RankCompressor.hpp"
#include "ApproxSVD.hpp"

namespace ACE{

template <typename T>
class RankCompressor_ApproxSVD_ScalarType: public RankCompressor_ScalarType<T>{
public:
  double threshold;
  double eps2;
  double sum_threshold;
  int maxk;
  bool precondition_repeat;
  bool forceQR;
  int count, DUMP_SVD;

  ApproxSVD<T> svd_bck;

  //increase threshold logarithmically from 'threshold' to 'threshold_to' within 'Nrange' SVDs
  double threshold_to;
  int Nrange;
 

  // check if parameters are set up so that any compression can happen
  bool has_effect()const;

  //choose how many singular values should be kept:  
  int get_new_dim(const Eigen::VectorXd &sv);
   
  void compress(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A, 
                         Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &L, 
                         Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &R, 
                         bool low_to_high);

  void setup(Parameters &param);

  inline RankCompressor_ApproxSVD_ScalarType(Parameters &param){
    setup(param);
  }
  inline RankCompressor_ApproxSVD_ScalarType(double thresh=0., double sum_thr_=0., int maxk_=0, int reortho=0, int dump=-1)
   : threshold(thresh),sum_threshold(sum_thr_), maxk(maxk_), count(0), DUMP_SVD(dump){
  }
  inline virtual ~RankCompressor_ApproxSVD_ScalarType(){
  }
};

typedef RankCompressor_ApproxSVD_ScalarType<std::complex<double> > RankCompressor_ApproxSVD;
typedef RankCompressor_ApproxSVD_ScalarType<double> RankCompressor_ApproxSVD_real;


}//namespace
#endif
