#ifndef ACE_SWEEP_TRAFO_PROCESSOR_DEFINED_H
#define ACE_SWEEP_TRAFO_PROCESSOR_DEFINED_H

#include "MPS_Matrix.hpp"

/* Abstraction to track the effects of SVD compression transformations, e.g.,
on environment operators */
namespace ACE{

template <typename T> 
class Sweep_Trafo_Processor_ScalarType{
public:

  virtual void process_low_to_high(int n, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &R);

  virtual void process_high_to_low(int n, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &L);

};

typedef Sweep_Trafo_Processor_ScalarType<std::complex<double> > Sweep_Trafo_Processor;
typedef Sweep_Trafo_Processor_ScalarType<double> Sweep_Trafo_Processor_real;

}//namespace
#endif
