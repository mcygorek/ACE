#include "RealFunction_Interpolate.hpp"
#include "Function.hpp"
#include "ReadTable.hpp"
#include <vector>
#include <algorithm>
#include <iostream>

namespace ACE{

  bool RealFunction_Interpolate::cmp_less(const std::pair<double,double> &p1,
                                          const std::pair<double,double> &p2){
    return p1.first<p2.first;
  }

  void RealFunction_Interpolate::sort(){
    std::sort(val.begin(), val.end(), cmp_less);
  }

  void RealFunction_Interpolate::read(const std::string &fname, int col1, int col2){
    val.clear();
    ReadTable tab(fname, col1, col2);
    for(size_t i=0; i<tab.size(); i++){
      val.push_back(std::make_pair(tab[i][0],tab[i][1]));
    }
    if(val.size()<1){
      std::cerr<<"Couldn't read a single line from file '"<<fname<<"'!"<<std::endl;  
      exit(1);
    }
    sort();
  }

  double RealFunction_Interpolate::f(double x)const{
    if(val.size()<1)return 0.;
    if(x<=val[0].first)return val[0].second;
    if(x>=val.back().first)return val.back().second;
    
    std::vector<std::pair<double,double> >::const_iterator it=
      std::lower_bound(val.begin(), val.end(), std::make_pair(x,0.), cmp_less);
    
    if(it==val.begin())return val[0].second;
    if(it==val.end())return val.back().second;

    std::vector<std::pair<double,double> >::const_iterator it2=it;
    it--;

    double dx=it2->first - it->first;
    if(fabs(dx)<1e-20){
      std::cerr<<"Cannot interpolate between sampling points with identical x values!"<<std::endl;
      exit(1);
    }
    double frac=(x - it->first)/dx;
    return it->second*(1.-frac)+it2->second*frac;
  }

}//namespace
