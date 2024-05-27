#ifndef RANK_COMPRESSOR_DEFINED_H
#define RANK_COMPRESSOR_DEFINED_H

#include "MPS.hpp"
#include <Eigen/Core>
#include "Sweep_Trafo_Processor.hpp"
#include "Compress_Trafo_At.hpp"

namespace ACE{
class Parameters;

template <typename T>
class RankCompressor_ScalarType{
public:
  int debug;
  bool keep_weight;

  //Note: left_to_right=true means: start from a[0]
  virtual void compress(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A, 
                        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &L, 
                        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &R, 
                        bool low_to_high)=0;


  virtual void compress(MPS_Matrix_ScalarType<T> &a, 
               Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &L, 
               Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &R, 
               bool low_to_high, double s=0.);
  
  virtual void compress_keep_largest(
                     MPS_Matrix_ScalarType<T> &a, 
                     Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &L, 
                     Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &R,  
                     bool low_to_high);

  void sweep_block_low_to_high(int n, MPS_ScalarType<T> &mps,
      double keep_weight=0., Sweep_Trafo_Processor_ScalarType<T> *proc=NULL);

  void sweep_block_high_to_low(int n, MPS_ScalarType<T> &mps,
      double keep_weight=0., Sweep_Trafo_Processor_ScalarType<T> *proc=NULL);

  void low_end_multiply_and_compress(MPS_ScalarType<T> &mps, 
                  const MPS_ScalarType<T> & other, int sweep_start=0);

  void high_end_multiply_and_compress(MPS_ScalarType<T> & mps,
            const MPS_ScalarType<T> & other, double keep_weight=0., 
                                        Compress_Trafo_At *cta=NULL);

  void multiply_and_compress(MPS_ScalarType<T> & mps, 
                     const MPS_ScalarType<T> & other);

  virtual void sweep(MPS_ScalarType<T> & mps, double keep_weight=0.);

  inline RankCompressor_ScalarType(){
    debug=0;
    keep_weight=true;
  }
  inline virtual ~RankCompressor_ScalarType(){
  }
};

typedef RankCompressor_ScalarType<std::complex<double> > RankCompressor;
typedef RankCompressor_ScalarType<double> RankCompressor_real;


template <typename T>
class RankCompressor_None_ScalarType: public RankCompressor_ScalarType<T>{
public:
  virtual void compress(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A, 
                        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &L, 
                        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &R, 
                        bool low_to_high);
};


typedef RankCompressor_None_ScalarType<std::complex<double> > RankCompressor_None;
typedef RankCompressor_None_ScalarType<double> RankCompressor_None_real;


}//namespace
#endif
