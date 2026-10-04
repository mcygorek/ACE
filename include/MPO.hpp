#ifndef MPO_DEFINED_H_
#define MPO_DEFINED_H_

#include "FourLeg.hpp"
#include "TruncationLayout.hpp"

namespace ACE{

template <typename T> class MPO_ScalarType{
public:
  std::vector<FourLeg_ScalarType<T> > a;
  
  inline FourLeg_ScalarType<T> & operator[](int l){
    if(l<0||l>=a.size()){
      std::cerr<<"VMPO: operator["<<l<<"] out of bounds (a.size()="<<a.size()<<")!"<<std::endl;
      exit(1);
    }
    return a[l];
  }
  inline const FourLeg_ScalarType<T> & operator[](int l)const{
    if(l<0||l>=a.size()){
      std::cerr<<"VMPO: operator["<<l<<"] out of bounds (a.size()="<<a.size()<<")!"<<std::endl;
      exit(1);
    }
    return a[l];
  }
  inline FourLeg_ScalarType<T> & back(){
    if(a.size()<1){
      std::cerr<<"VMPO: operator["<<0<<"] out of bounds (a.size()="<<a.size()<<")!"<<std::endl;
      exit(1);
    }
    return a[a.size()-1];
  }
  inline const FourLeg_ScalarType<T> & back()const{
    if(a.size()<1){
      std::cerr<<"VMPO: operator["<<0<<"] out of bounds (a.size()="<<a.size()<<")!"<<std::endl;
      exit(1);
    }
    return a[a.size()-1];
  }

  inline void scale_first(const T & scale){ 
    if(a.size()>0){ a[0].scale(scale); }
  }

  inline void clear(){ a.clear(); }

  size_t size()const{ return a.size(); }

  void resize(int n, const FourLeg_ScalarType<T> & obj=FourLeg_ScalarType<T>() ){ a.resize(n, obj); }

  void sweep_single_d2_d4(const TruncatedSVD &trunc, int l, int verbosity=1);
  void sweep_d2_d4(const TruncatedSVD &trunc, int verbosity=1);

  void sweep_single_d4_d2(const TruncatedSVD &trunc, int l, int verbosity=1);
  void sweep_d4_d2(const TruncatedSVD &trunc, int verbosity=1);

  static MPO_ScalarType<T> Identity_d1_d3(const std::vector<int> & dims);

  void print_dims(std::ostream &os=std::cout)const;

  inline void copy(const MPO_ScalarType<T> &other){
    a.clear();
    a.resize(other.a.size());
    for(size_t i=0; i<a.size(); i++){
      a[i].copy(other.a[i]);
    }
  }
  inline void swap(MPO_ScalarType<T> &other){ a.swap(other.a); }

  inline MPO_ScalarType &operator=(const MPO_ScalarType<T> &other){ copy(other); return *this; }
  inline MPO_ScalarType(const MPO_ScalarType<T> &other){ copy(other); }
  MPO_ScalarType(){}
  MPO_ScalarType(int n, const FourLeg_ScalarType<T> & obj=FourLeg_ScalarType<T>() ){resize(n, obj);}
  ~MPO_ScalarType(){}
};
typedef MPO_ScalarType<std::complex<double> > MPO;
typedef MPO_ScalarType<double> MPO_real;


}//namespace
#endif 
