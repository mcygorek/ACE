#pragma once
#ifndef ACE_READ_TEMPERATURE_DEFINED_H_
#define ACE_READ_TEMPERATURE_DEFINED_H_

#include "Parameters.hpp"
#include "Constants.hpp"

namespace ACE{

extern double readTemperature(Parameters &param, const std::string &prefix="");

}//namespace
#endif
