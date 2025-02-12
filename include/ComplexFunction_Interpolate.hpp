#ifndef ACE_COMPLEXFUNCTION_INTERPOLATE_DEFINED_H
#define ACE_COMPLEXFUNCTION_INTERPOLATE_DEFINED_H

#include "Function.hpp"
#include "ReadTable.hpp"
#include <vector>
#include <algorithm>
#include <iostream>

namespace ACE{

class ComplexFunction_Interpolate: public ComplexFunction{
public:
  std::vector<std::pair<double,std::complex<double> > > val; // (x_i, f(x_i)); assume sorted

  static bool cmp_less(const std::pair<double,std::complex<double> > &p1,
                       const std::pair<double,std::complex<double> > &p2);

  inline void sort(){
    std::sort(val.begin(), val.end(), cmp_less);
  }
  void read(const std::string &fname, int col1=0, int col2=1, int col3=2);

  virtual std::complex<double> f(double x)const;

  inline ComplexFunction_Interpolate(const std::vector<std::pair<double,std::complex<double> > > &val_) : val(val_){
  }
  
  inline ComplexFunction_Interpolate(const std::string &fname, int col1=0, int col2=1, int col3=2){
    read(fname, col1, col2, col3);
  }
  inline ComplexFunction_Interpolate(){}
};

}//namespace

#endif
