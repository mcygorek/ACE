#include "ACE.hpp"
#include "Parameters.hpp"
#include "ReadTable.hpp"
#include <Eigen/Dense>
#include "LeastSquares.hpp"

using namespace ACE;

Eigen::VectorXd f(Eigen::VectorXd v){ //exp(-k*x) ; v(0)=k
  Eigen::VectorXd ret(10);
  for(int i=0; i<ret.rows(); i++){
    ret(i)=exp(-v(0)*i);
  }
  return ret;
}
Eigen::MatrixXd J(Eigen::VectorXd v){
  Eigen::MatrixXd ret(10,1);
  for(int i=0; i<ret.rows(); i++){
    ret(i,0)=-((double)i)*exp(-v(0)*i);
  }
  return ret;
}

int main(int args, char** argv){

  Parameters param(args, argv);
//  ReadTable

  Eigen::VectorXd guess(1); guess(0)=2.;
  Eigen::VectorXd ref=f(guess);
  guess(0)=1;
  
  guess=LeastSquares::GaussNewton(ref, f, J, guess, 1e-8, 100, 1);
/*
  for(int loop=0; loop<5; loop++){
    std::cout<<"loop="<<loop<<std::endl;
    std::pair<Eigen::VectorXd,Eigen::VectorXd> retGS =
                             LeastSquares::GaussNewton_step(ref, f, J, guess);
    guess=retGS.first;
    std::cout<<"guess="<<guess.transpose()<<std::endl;
    std::cout<<"residual="<<retGS.second.norm()<<std::endl;
  }
*/

  std::cout<<"Final: "<<guess.transpose()<<std::endl;

  return 0;
}


