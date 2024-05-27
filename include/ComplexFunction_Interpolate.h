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
                       const std::pair<double,std::complex<double> > &p2){
    return p1.first<p2.first;
  }

  void sort(){
    std::sort(val.begin(), val.end(), cmp_less);
  }
  void read(const std::string &fname, int col1=0, int col2=1, int col3=2){
    val.clear();
    ReadTable tab(fname, col1, col2, col3);
    for(size_t i=0; i<tab.size(); i++){
//      std::cout<<tab[i][0]<<" "<<std::complex<double>(tab[i][1],tab[i][2])<<std::endl;
      val.push_back(std::make_pair(tab[i][0],
        std::complex<double>(tab[i][1],tab[i][2])));
    }
    sort();
  }

  virtual std::complex<double> f(double x)const{
    if(val.size()<1)return 0.;
    if(x<=val[0].first)return val[0].second;
    if(x>=val.back().first)return val.back().second;
    
    std::vector<std::pair<double,std::complex<double> > >::const_iterator it=
      std::lower_bound(val.begin(), val.end(), std::make_pair(x,0.), cmp_less);
    
    if(it==val.begin())return val[0].second;
    if(it==val.end())return val.back().second;

    std::vector<std::pair<double,std::complex<double> > >::const_iterator it2=it;
    it--;

    double dx=it2->first - it->first;
    if(fabs(dx)<1e-20){
      std::cerr<<"Cannot interpolate between sampling points with identical x values!"<<std::endl;
      exit(1);
    }
    double frac=(x - it->first)/dx;
    return it->second*(1.-frac)+it2->second*frac;
  }

  ComplexFunction_Interpolate(const std::string &fname, int col1=0, int col2=1, int col3=2){
    read(fname, col1, col2, col3);
  }
  ComplexFunction_Interpolate(){}
};

}//namespace

#endif
