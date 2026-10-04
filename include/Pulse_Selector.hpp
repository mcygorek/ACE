#ifndef ACE_PULSE_SELECTOR_DEFINED_H
#define ACE_PULSE_SELECTOR_DEFINED_H

#include "Function.hpp"

namespace ACE{

std::pair<ComplexFunctionPtr, Eigen::MatrixXcd>  Pulse_Selector(
                                      const std::vector<std::string> &toks);

}//namespace
#endif
