#include "EDM_Index.hpp"
#include "DummyException.hpp"
#include <iostream>

namespace ACE {

bool EDM_Index::operator<(const EDM_Index &other)const{
  if(other.size() != size()){
    std::cerr<<"EDM_Index::operator<: Comparing index lists of different lengths!"<<std::endl;
    throw DummyException();
  }
  for(int i=0; i<list.size(); i++){
    if(list[i]<other[i])return true;
    if(list[i]>other[i])return false;
  }
  return false;
}
bool EDM_Index::operator==(const EDM_Index &other)const{
  if(other.size() != size()){
    return false;
//    std::cerr<<"EDM_Index::operator==: Comparing index lists of different lengths!"<<std::endl;
//    throw DummyException();
  }
  for(int i=0; i<list.size(); i++){
    if(list[i]!=other[i])return false;
  }
  return true;
}



void EDM_Index::increase(const EDM_Index & Ldim){
  check_rank(Ldim.size());
  for(int i=(int)size()-1; i>=0; i--){
    list[i]++;
    if(list[i]<Ldim[i]){
      break;
    }else{
      if(i!=0){
        list[i]=0;
      }
    }
  }
}
void EDM_Index::increase_at(int i, const EDM_Index & Ldim){
  check_rank(Ldim.size());
  set_zero_from(i+1);
  for(int j=i; j>=0; j--){
    list[j]++;
    if(list[j]<Ldim[j]){
      break;
    }else{
      if(j!=0){
        list[j]=0;
      }
    }
  }
}

void EDM_Index::print(std::ostream &os)const{
  for(int i=0; i<list.size(); i++){ 
    if(i>0)os<<" ";
    os<<list[i];
  }
}

void EDM_Index::check_rank(int rank)const{
  if(size()!=rank){
    std::cerr<<"EDM_Index: Mismatching ranks "<<size()<<" vs. "<<rank<<"!"<<std::endl;
    throw DummyException();
  }
}
void EDM_Index::check_in_range(int site)const{
  if(site<0||site>=list.size()){
    std::cerr<<"EDM_Index: System index out of range: "<<site<<" vs. [0:"<<list.size()<<"["<<std::endl;
    throw DummyException();
  }
}
void EDM_Index::check_in_range(int site, int index)const{
  check_in_range(site);
  if(index<0||index>=list[site]){
    std::cerr<<"EDM_Index: index at site "<<site<<" out of range: "<<index<<" vs. [0:"<<list[site]<<"["<<std::endl;
    throw DummyException();
  }
}
void EDM_Index::check_in_range(const EDM_Index &I)const{
  if(list.size()!=I.size()){
    std::cerr<<"EDM_Index::check_in_range(I): list.size()!=I.size()"<<std::endl;
    throw DummyException();
  }
  for(int r=0; r<I.size(); r++){
    if(I[r]<0||I[r]>=list[r]){
      std::cerr<<"EDM_Index::check_in_range(Ldim): I["<<r<<"] out of bounds"<<std::endl;
      throw DummyException();
    }
  }  
}

void EDM_Index::check_if_equal(const EDM_Index &I)const{
  if(list.size()!=I.size()){
    std::cerr<<"EDM_Index::check_if_equal(I): list.size()!=I.size()"<<std::endl;
    throw DummyException();
  }
  for(int r=0; r<I.size(); r++){
    if(I[r]!=list[r]){
      std::cerr<<"EDM_Index::check_of_equal: I["<<r<<"] != list["<<r<<"]"<<std::endl;
      throw DummyException();
    }
  }  
}
}//namespace


std::size_t std::hash<ACE::EDM_Index>::operator()(const ACE::EDM_Index& I)const noexcept{
  std::size_t seed=0;
  for(int i=0; i<I.list.size(); i++){
    int i_hash=std::hash<int>{}(I.list[i]);
    seed ^= i_hash + 0x9e3779b9 + (seed<<6) + (seed>>2);
    }
  return seed;
}



