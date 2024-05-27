#pragma once
#ifndef FT_PARAMETERS_DEFINED_H
#define FT_PARAMETERS_DEFINED_H

#include "slowFT.hpp"

namespace ACE{

struct FT_Parameters{
  double wa, we, dt, ta;
  int Nsubdiv, Ndiscr;
  int sig;
  
  FT_Parameters(){
    sig=-1; Nsubdiv=10;
  }
};

inline std::vector<std::complex<double> > slowFT(const std::vector<std::complex<double> > & in, const FT_Parameters & ftp){
  return slowFT(in, ftp.wa, ftp.we, ftp.Ndiscr, ftp.ta, ftp.dt, ftp.sig, ftp.Nsubdiv);
}

}//namespace
#endif
