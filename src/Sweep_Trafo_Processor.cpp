#include "Sweep_Trafo_Processor.hpp"
#include "MPS_Matrix.hpp"

/* Abstraction to track the effects of SVD compression transformations, e.g.,
on environment operators */
namespace ACE{

template <typename T> void Sweep_Trafo_Processor_ScalarType<T>::
    process_low_to_high(int n, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &R){};

template <typename T> void Sweep_Trafo_Processor_ScalarType<T>::
    process_high_to_low(int n, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &L){};


template class Sweep_Trafo_Processor_ScalarType<std::complex<double> >;
template class Sweep_Trafo_Processor_ScalarType<double>;

}//namespace
