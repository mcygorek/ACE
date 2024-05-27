#ifndef N_CHOOSE_K_DEFINED_H
#define N_CHOOSE_K_DEFINED_H

#include <vector>

namespace ACE{
//Calculate n choose k efficiently by precomputing values (lazy evaluation)
class N_Choose_K{
  std::vector<std::vector<int> > precomp;
public:

  void precompute(int n);

  int operator()(int n, int k);
  
  inline N_Choose_K(int nmax=0){
    if(nmax>0)precompute(nmax);
  }

};
}//namespace
#endif
