#ifndef ACE_ABSTRACT_EIGEN_MPS_DEFINED_H
#define ACE_ABSTRACT_EIGEN_MPS_DEFINED_H

#include "EigenMPS_Matrix.h"

class Abstract_EigenMPS{
public: 
  virtual int size()const=0;
   
  virtual const EigenMPS_Matrix & operator[](int n)const=0;

  void print_dims(std::ostream &os=std::cout) const{
    for(size_t n=0; n<size(); n++){
      operator[](n).print_dims(os);
      os<<std::endl;
    }
  }
  void print_dims(const std::string &filename) const{
    std::ofstream ofs(filename.c_str());
    print_dims(ofs);
  }
  int get_max_dim()const{
    int max_dim=0;
    for(size_t n=0; n<size(); n++){
      if(operator[](n).get_dim_d2()>max_dim)max_dim=operator[](n).get_dim_d2();
    }
    return max_dim;
  }
  int get_max_dim_pos()const{
    int max_dim=0;
    int pos=-1;
    for(size_t n=0; n<size(); n++){
      if(operator[](n).get_dim_d2()>max_dim){
        max_dim=operator[](n).get_dim_d2();
        pos=n;
      }
    }
    return pos;
  }
  void print_max_dim(std::ostream &os=std::cout) const{
    os<<get_max_dim();
  }
  virtual void print_inner_dims(std::ostream &os=std::cout) const{
    for(int n=0; n<size(); n++){
      if(n==0){
        os<<operator[](n).get_dim_d1()<<" "<<operator[](n).get_dim_d2();
      }else{
        os<<" "<<operator[](n).get_dim_d2();
      }
    }
  }
  double max_element_abs()const{
    double max=0;
    for(size_t n=0; n<size(); n++){
      double m=operator[](n).max_element_abs();
      if(m>max)max=m;
    }
    return max;
  }


  virtual ~Abstract_EigenMPS(){}
};

#endif
