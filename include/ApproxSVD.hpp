#ifndef ACE_APPROX_SVD_DEFINED_H
#define ACE_APPROX_SVD_DEFINED_H

//#include "Eigen_fwd.hpp"
#include <Eigen/Dense>
#include <vector>

namespace ACE{

template<typename ScalarType> double conditionalReal(ScalarType z);
template<> double conditionalReal<std::complex<double> >(std::complex<double> z);


template<typename ScalarType> class ApproxSVD{
public:
  typedef Eigen::Matrix<ScalarType, -1, -1, 0, -1, -1> MatrixType;
  MatrixType L,M,R;  //left, middle and right matrices. 
  MatrixType QRremainder; 
 
  
  void rotate_single(int p, int q);

  void initialize(const MatrixType &A, bool forceQR=false);

  MatrixType reconstructL();
  MatrixType reconstructR();

  static bool sort_second_descending(const std::pair<int, double> & first, 
                                     const std::pair<int, double> & second);
  
  void sort();
  
  void calculate(double eps, int max_loop=-1);
  
  Eigen::VectorXd singularValues()const;
  
  void swap(ApproxSVD &other);

  void compute(const MatrixType &A, double eps, bool forceQR=false);
  
  inline ApproxSVD(){}
  inline ApproxSVD(const MatrixType &A, double eps, bool forceQR=false){
    compute(A, eps, forceQR);
  }
  inline virtual ~ApproxSVD(){}
};


template<> void ApproxSVD<std::complex<double> >::rotate_single(int p, int q);


}//namespace
#endif
