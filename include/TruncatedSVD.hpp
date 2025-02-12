#ifndef ACE_TRUNCATED_SVD_DEFINED_H
#define ACE_TRUNCATED_SVD_DEFINED_H

#include "Eigen_fwd.hpp"
#include "Parameters.hpp"
#include "SelectIndices.hpp"
#include "PassOn.hpp"
#include <string>

namespace ACE{


//TODO: need more general structure that can also deal with (RR)QR:
//Ideally allowing to extract a "pass_on", with the inverse realizable via 
//a .solve(..) call, as well as a "keep_here"
//Also needs more general notion of "weights" rather then the SVDs
template <typename T> struct TruncatedSVD_Return_T{
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> U;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Vdagger;
  Eigen::VectorXd sigma;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Residual;
};
typedef TruncatedSVD_Return_T<std::complex<double> > TruncatedSVD_Return;
typedef TruncatedSVD_Return_T<double> TruncatedSVD_Return_real;


template <typename T> class TruncatedSVD_T{
public:

  //relevant for concrete sweeps:
  double threshold;
  int    maxk;
  double keep;  //scaling factor: keep <= 0: keep largest singular value.

  bool use_QR; //Orthogonalization using QR instead of compressing sweep
  double Tikhonov; //Tikhonov regularization factor for inverse matrix
  /*
  To obtain a concrete threshold/maxk, we consider N lines, with M sweeps each
  
  */

  virtual int get_truncated_dim(const Eigen::VectorXd &svals)const;
  SelectIndices get_select_indices(const Eigen::VectorXd & s1, const Eigen::VectorXd & s2)const;

//Perform truncated SVD (returns U, sigma, Vdagger):
  virtual TruncatedSVD_Return_T<T> compress(
       const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & A,
       bool calculate_residual=false) const;

//Compression for forward sweep (updates argument, returns PassOn)
//[is abstract enough to also process QR factorizations]:
  virtual PassOn_T<T> compress_forward(
             Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & A, 
             Eigen::VectorXd & weights) const;
//same for backward sweep:
  virtual PassOn_T<T> compress_backward(
             Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & A, 
             Eigen::VectorXd & weights) const;
 

  virtual void scale(double threshold_ratio, double maxk_ratio);

  virtual bool do_compress()const;

  virtual void print_info(std::ostream &ofs=std::cout)const;
  virtual void setup(Parameters &param);

  TruncatedSVD_T(Parameters &param){
    setup(param);
  }
  TruncatedSVD_T(double threshold_=0, int maxk_=0, int keep_=0, bool use_QR_=false, double Tikhonov_=0.){
    threshold=threshold_;
    maxk=maxk_;
    keep=keep_;
    use_QR=use_QR_;
    Tikhonov=Tikhonov_;
  }
};


typedef TruncatedSVD_T<std::complex<double> > TruncatedSVD;
typedef TruncatedSVD_T<double> TruncatedSVD_real;

}//namespace
#endif
