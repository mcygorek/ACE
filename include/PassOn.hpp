#ifndef ACE_PASS_ON_DEFINED_H
#define ACE_PASS_ON_DEFINED_H
#include <Eigen/Dense>

//structure to contain matrices to be passed on during sweeps:
template <typename T> struct PassOn_T{
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> P; //R/L
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Pinv; //Rinv/Linv

  inline void set(int d){
    Pinv=P=Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(d, d);
  }
  inline PassOn_T<T> get_reversed()const{
    PassOn_T<T> p;   
    p.P=Pinv;
    p.Pinv=P;
    return p;
  }
  inline PassOn_T(int d=0){
    set(d);
  }
};

typedef PassOn_T<std::complex<double> > PassOn;
typedef PassOn_T<double> PassOn_real;

#endif
