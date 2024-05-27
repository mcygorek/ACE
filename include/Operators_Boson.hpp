#pragma once
#ifndef OPERATORS_BOSON_DEFINED_H
#define OPERATORS_BOSON_DEFINED_H

#include <Eigen/Core>
#include "Operators.hpp"

namespace ACE{
namespace Operators_Boson{

extern Eigen::MatrixXcd a(int n_max);
extern Eigen::MatrixXcd adagger(int n_max);
extern Eigen::MatrixXcd n(int n_max);
extern Eigen::MatrixXcd id(int N);
extern Eigen::MatrixXcd vacuum(int n_max);

};

namespace Operators_Boson_Offset{

extern Eigen::MatrixXcd a(int base, int M);
extern Eigen::MatrixXcd adagger(int base, int M);
extern Eigen::MatrixXcd n(int base, int M);
extern Eigen::MatrixXcd lowest(int M);
};

}//namespace
#endif
