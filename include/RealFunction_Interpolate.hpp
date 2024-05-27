#ifndef ACE_REALFUNCTION_INTERPOLATE_DEFINED_H
#define ACE_REALFUNCTION_INTERPOLATE_DEFINED_H

#include "Function.hpp"
#include "ReadTable.hpp"
#include <vector>
#include <algorithm>

namespace ACE{

class RealFunction_Interpolate: public RealFunction{
public:
  std::vector<std::pair<double,double> > val; // (x_i, f(x_i)); assume sorted

  static bool cmp_less(const std::pair<double,double> &p1,
                       const std::pair<double,double> &p2);
  void sort();

  void read(const std::string &fname, int col1=0, int col2=1);

  virtual double f(double x)const;

  inline RealFunction_Interpolate(const std::string &fname, int col1=0, int col2=1){
    read(fname, col1, col2);
  }
  inline RealFunction_Interpolate(){}
};


}//namespace
#endif
