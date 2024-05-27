#include "LeastSquares.hpp"
#include <Eigen/SVD>
#include <iostream>
#include <cstdlib>


namespace ACE{
namespace LeastSquares{

std::pair<Eigen::VectorXd, Eigen::VectorXd>  GaussNewton_step(
    Eigen::VectorXd reference,
    std::function<Eigen::VectorXd(Eigen::VectorXd)> function, 
    std::function<Eigen::MatrixXd(Eigen::VectorXd)> jacobian,
    Eigen::VectorXd guess){

    //for(int loop=0; loop<5; loop++)
    Eigen::VectorXd f = function(guess);
    if(f.rows()!=reference.rows()){
      std::cerr<<"GaussNewton: f.rows()!=reference.rows()!"<<std::endl;
      exit(1);
    }
    Eigen::VectorXd r = reference-f;
    Eigen::MatrixXd J = jacobian(guess);
    if(J.cols()!=guess.rows()){
      std::cerr<<"GaussNewton: J.cols()="<<J.cols()<<"!=guess.rows()="<<guess.rows()<<"!"<<std::endl;
      exit(1);
    }
    if(J.rows()!=r.rows()){
      std::cerr<<"GaussNewton: J.rows()!=r.rows()!"<<std::endl;
      exit(1);
    }
    Eigen::MatrixXd A=J.transpose()*J;

//    Eigen::JacobiSVD<Eigen::MatrixXd>
//                    svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
//    Eigen::VectorXd sinv=svd.singularValues();
//    std::cout<<"singular values: "<<sinv.transpose()<<std::endl;
//    for(int i=0; i<sinv.rows(); i++){
//      sinv(i)=1./sinv(i);
//    }
//    Eigen::MatrixXd Ainv=svd.matrixV()*(sinv.asDiagonal())*svd.matrixU().adjoint();     
//    guess+=Ainv*J.transpose()*r;

/*
  std::cout<<"test0"<<std::endl;
  std::cout<<"guess.transpose()="<<guess.transpose()<<std::endl;
  std::cout<<"test1"<<std::endl;
  std::cout<<"A.col(0).transpose()="<<A.col(0).transpose()<<std::endl;
  std::cout<<"test2"<<std::endl;
  std::cout<<"J.row(0)="<<J.row(0)<<std::endl;
  std::cout<<"test2.5"<<std::endl;
  std::cout<<"r.head(10).transpose()="<<r.head(10).transpose()<<std::endl;
  std::cout<<"test3"<<std::endl;
  std::cout<<"J.cols()="<<J.cols()<<" J.rows()="<<J.rows()<<" r.rows()="<<r.rows()<<std::endl;
  Eigen::VectorXd Jtr=J.transpose()*r;
  std::cout<<"Jtr.transpose()="<<Jtr.transpose()<<std::endl;
  std::cout<<"test3.5"<<std::endl;
*/
  guess+=A.colPivHouseholderQr().solve(J.transpose()*r);
//  std::cout<<"guess.transpose()="<<guess.transpose()<<std::endl;
//  std::cout<<"test4"<<std::endl;

  return std::make_pair(guess, r);
}


Eigen::VectorXd GaussNewton(
    Eigen::VectorXd reference,
    std::function<Eigen::VectorXd(Eigen::VectorXd)> function,
    std::function<Eigen::MatrixXd(Eigen::VectorXd)> jacobian,
    Eigen::VectorXd guess,
    double epsilon, int maxloop, int verbose){

  double last_residual=0;
  for(int loop=0; loop<maxloop; loop++){
    if(verbose>0)std::cout<<"loop="<<loop<<std::endl;
    std::pair<Eigen::VectorXd,Eigen::VectorXd> retGS =
           LeastSquares::GaussNewton_step(reference, function, jacobian, guess);
    guess=retGS.first;
    if(verbose>0)std::cout<<"guess="<<guess.transpose()<<std::endl;
    double residual=retGS.second.norm();
    if(verbose>0)std::cout<<"residual="<<residual<<std::endl;

    if(loop>0 && fabs(last_residual-residual)<epsilon){break;}
    last_residual=residual;
  }
  return guess;
}

}//namespace
}//namespace
