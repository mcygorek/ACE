#ifndef ACE_TRUNCATION_RETURN_DEFINED_H
#define ACE_TRUNCATION_RETURN_DEFINED_H


namespace ACE{

template <typename T> struct TruncationReturn_T{
  PassOn pass_On;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> keep;
  Eigen::VectorXd weights;
};
typedef TruncationReturn_T<std::complex<double> > TruncationReturn;
typedef TruncationReturn_T<double> TruncationReturn_real;



}
#endif
