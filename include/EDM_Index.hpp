#ifndef ACE_EDM_INDEX_DEFINED_H
#define ACE_EDM_INDEX_DEFINED_H

#include <vector>
#include <cstdlib>
#include <iostream>

namespace ACE {

class EDM_Index{
public:
  std::vector<int> list;

  inline void clear(){list.clear();}
  inline void resize(size_t r){list.resize(r);}
  inline void erase(int r){list.erase(list.begin()+r);}
  inline size_t size()const{return list.size();}
  inline int & operator[](int i){return list[i];}
  inline const int & operator[](int i)const{return list[i];}
  inline void swap(EDM_Index & other){
    list.swap(other.list);
  }
  bool operator<(const EDM_Index &other)const;
  bool operator==(const EDM_Index &other)const;

  //looping
  void increase(const EDM_Index & Ldim);
  //increases index i by one and sets indices after that to zero
  void increase_at(int i, const EDM_Index & Ldim);

  void print(std::ostream &os=std::cout)const;

  //checks assuming *this is Ldim
  void check_rank(int rank)const;
  void check_in_range(int site)const;
  void check_in_range(int site, int index)const;
  void check_in_range(const EDM_Index &I)const;
  void check_if_equal(const EDM_Index &other)const;
  

  inline void set_zero(int rank){
    list=std::vector<int>(rank,0);
  }
  inline void set_zero_from(int i){
    for(int j=i; j<size(); j++){
      list[j]=0;
    }
  }

  EDM_Index(){}
  EDM_Index(int rank){
    set_zero(rank);
  }
};
}

template<> struct std::hash<ACE::EDM_Index>{
    std::size_t operator()(const ACE::EDM_Index& I) const noexcept;
};


#endif
