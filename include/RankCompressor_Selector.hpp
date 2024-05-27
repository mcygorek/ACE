#ifndef RANK_COMPRESSOR_SELECTOR_DEFINED_H
#define RANK_COMPRESSOR_SELECTOR_DEFINED_H

#include "Smart_Ptr.h"
#include "RankCompressor.hpp"

namespace ACE{
class Parameters;

typedef Smart_Ptr<RankCompressor> RankCompressor_Ptr;

RankCompressor_Ptr RankCompressor_Selector(Parameters &param, bool verbose=false);

}//namespace
#endif
