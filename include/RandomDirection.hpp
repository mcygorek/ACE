#ifndef RANDOM_DIRECTION_DEFINED_H
#define RANDOM_DIRECTION_DEFINED_H

#include "Eigen_fwd.hpp"

/** Get random direction to initialize spin baths */
namespace ACE{

double random_double();

Eigen::Vector3d RandomDirection();

void RandomAngles(double &theta, double &phi);

}//namespace
#endif
