#ifndef ACE_RANKCOMPRESSOR_RRQR_DEFINED_H
#define ACE_RANKCOMPRESSOR_RRQR_DEFINED_H

#include "RankCompressor.hpp"

namespace ACE{

class RankCompressor_RRQR: public RankCompressor{
public: 
  double threshold;
  
  virtual void compress(Eigen::MatrixXcd &A, Eigen::MatrixXcd &L, Eigen::MatrixXcd &R, bool low_to_high);

  inline RankCompressor_RRQR(double thresh): threshold(thresh){
  }
  inline virtual ~RankCompressor_RRQR(){
  }
};

}//namespace
#endif
