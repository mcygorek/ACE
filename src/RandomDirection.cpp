#include "RandomDirection.hpp"
#include <cmath>
#include <Eigen/Core>

/** Get random direction to initialize spin baths */
namespace ACE{

double random_double(){
  Eigen::MatrixXd m=Eigen::MatrixXd::Random(1,1); //range: [-1:1]
  return (m(0,0)+1.)/2.; //range: [0:1]
}

Eigen::Vector3d RandomDirection(){
  Eigen::Vector3d v=Eigen::MatrixXd::Random(3,1);
  double n=v.squaredNorm();

  while(n<1e-6 || n>1.){
    v=Eigen::MatrixXd::Random(3,1);
    n=v.squaredNorm();
  } 
  return v/sqrt(n);
}
/*
Eigen::Vector3d RandomDirection(){
  Eigen::Vector3d v;
  double phi=2.*M_PI*random_double();
  double theta;
  do{
    theta=M_PI*random_double();
    v(0)=sin(theta)*cos(phi);
    v(1)=sin(theta)*sin(phi);
    v(2)=cos(theta);
  }while(random_double()>sin(theta));
  return v;
}
*/

void RandomAngles(double &theta, double &phi){
  Eigen::Vector3d v=RandomDirection();
  if(v(2)>=1.)v(2)=0.9999999999;
  theta=acos(v(2));
  phi=atan2(v(1), v(0));
}

}//namespace
