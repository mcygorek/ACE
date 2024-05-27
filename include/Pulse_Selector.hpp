#ifndef ACE_PULSE_SELECTOR_DEFINED_H
#define ACE_PULSE_SELECTOR_DEFINED_H

#include "Eigen_fwd.hpp"
#include "Function.hpp"
#include <vector>
#include <iosfwd>

namespace ACE{

std::pair<ComplexFunctionPtr, Eigen::MatrixXcd>  Pulse_Selector(
                                      const std::vector<std::string> &toks);

}//namespace
#endif
