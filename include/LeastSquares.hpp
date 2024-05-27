#ifndef ACE_LEAST_SQUARES_DEFINED_H
#define ACE_LEAST_SQUARES_DEFINED_H

#include <Eigen/Dense>
#include <functional>

namespace ACE{
namespace LeastSquares{

//single step of Gauss-Newton. Returns updated guess and residual
std::pair<Eigen::VectorXd, Eigen::VectorXd> GaussNewton_step(
    Eigen::VectorXd reference,
    std::function<Eigen::VectorXd(Eigen::VectorXd)> function,
    std::function<Eigen::MatrixXd(Eigen::VectorXd)> jacobian,
    Eigen::VectorXd guess);

Eigen::VectorXd GaussNewton(
    Eigen::VectorXd reference,
    std::function<Eigen::VectorXd(Eigen::VectorXd)> function,
    std::function<Eigen::MatrixXd(Eigen::VectorXd)> jacobian,
    Eigen::VectorXd guess,
    double epsilon, int maxloop, int verbose);


}//namespace
}//namespace
#endif
