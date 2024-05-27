#ifndef RANK_COMPRESSOR_POWERURV_DEFINED_H
#define RANK_COMPRESSOR_POWERURV_DEFINED_H

//#include "PowerURV.hpp"
#include "RankCompressor.hpp"
#include "Eigen_fwd.hpp"

namespace ACE{

class RankCompressor_PowerURV: public RankCompressor{
public: 
  double threshold;
  
  virtual void compress(Eigen::MatrixXcd &A, Eigen::MatrixXcd &L, Eigen::MatrixXcd &R, bool left_to_right);

  inline RankCompressor_PowerURV(double thresh): threshold(thresh){
  }
  inline virtual ~RankCompressor_PowerURV(){
  }
};

}//namespace

#endif
