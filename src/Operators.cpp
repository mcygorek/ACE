#include "PCH.hpp"
#include "Operators.hpp"

namespace ACE{

  Eigen::MatrixXcd sigma_x(){
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(2,2);
    mat(0,1)=1.;
    mat(1,0)=1.;
    return mat;
  }
  Eigen::MatrixXcd sigma_y(){
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(2,2);
    mat(0,1)=std::complex<double>(0.,1.);
    mat(1,0)=std::complex<double>(0.,-1.);
    return mat;
  }
  Eigen::MatrixXcd sigma_z(){
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(2,2);
    mat(0,0)=-1.;
    mat(1,1)=1.;
    return mat;
  }
  Eigen::MatrixXcd sigma_plus(){
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(2,2);
    mat(1,0)=1;
    return mat;
  }
  Eigen::MatrixXcd sigma_minus(){
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(2,2);
    mat(0,1)=1;
    return mat;
  }
}

