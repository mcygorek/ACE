#include "RankCompressor_RRQR.hpp"
#include "RRQR.hpp"
#include <Eigen/Dense>

namespace ACE{

  void RankCompressor_RRQR::compress(Eigen::MatrixXcd &A, Eigen::MatrixXcd &L, Eigen::MatrixXcd &R, bool low_to_high){
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

}//namespace
