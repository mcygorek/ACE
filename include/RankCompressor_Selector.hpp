#ifndef RANK_COMPRESSOR_SELECTOR_DEFINED_H
#define RANK_COMPRESSOR_SELECTOR_DEFINED_H

#include  <memory>
#include "RankCompressor.hpp"

namespace ACE{
class Parameters;

typedef std::shared_ptr<RankCompressor> RankCompressor_Ptr;

RankCompressor_Ptr RankCompressor_Selector(Parameters &param, bool verbose=false);

}//namespace
#endif
