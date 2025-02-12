#include "ComplexFunction_Interpolate.hpp"

namespace ACE{

bool ComplexFunction_Interpolate::cmp_less(
                       const std::pair<double,std::complex<double> > &p1,
                       const std::pair<double,std::complex<double> > &p2){
  return p1.first<p2.first;
}

void ComplexFunction_Interpolate::read(const std::string &fname, int col1, int col2, int col3){
  val.clear();
  ReadTable tab(fname, col1, col2, col3);
  for(size_t i=0; i<tab.size(); i++){
    val.push_back(std::make_pair(tab[i][0],
      std::complex<double>(tab[i][1],tab[i][2])));
  }
  sort();
}

std::complex<double> ComplexFunction_Interpolate::f(double x)const{
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


}//namespace

