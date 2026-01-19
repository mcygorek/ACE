#ifndef ACE_RANDOMIZED_COMPRESSION_DEFINED
#define ACE_RANDOMIZED_COMPRESSION_DEFINED

#include "ProcessTensorBuffer.hpp"

namespace ACE{

void Randomized_Combine(ProcessTensorBuffer & PTB, \
                        ProcessTensorBuffer & PTB2, \
                        int chi_new, const TruncatedSVD &trunc);

//void Randomized_test(ProcessTensorBuffer & PTB, \
                        ProcessTensorBuffer & PTB2, \
                        int chi_new, const TruncatedSVD &trunc);


}//namespace
#endif
