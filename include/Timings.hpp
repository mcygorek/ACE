#ifndef ACE_TIMINGS_DEFINED_H
#define ACE_TIMINGS_DEFINED_H
#include <chrono>

namespace ACE{

typedef std::chrono::time_point<std::chrono::steady_clock> time_point;

inline double time_diff(std::chrono::nanoseconds diff){
  return std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
}

inline time_point now(){
  return std::chrono::steady_clock::now();
}

}//namespace
#endif
