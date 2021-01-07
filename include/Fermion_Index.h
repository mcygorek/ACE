#ifndef FERMION_INDEX_DEFINED_H
#define FERMION_INDEX_DEFINED_H

#include <vector>
#include <fstream>
#include "N_Choose_K.h"

class Fermion_Index{
public:

  std::vector<int> list;
  int sig;

  int operator[](int i)const{return list[i];}
  int &operator[](int i){return list[i];}
  size_t size()const{return list.size();}
  void resize(int i){list.resize(i);}
  void push_back(int i){list.push_back(i);}
  void clear(){list.clear();}
  operator bool() const{return sig!=0;}
  Fermion_Index &operator<<(int i){
    list.push_back(i);
    return *this;
  }

  int sort(){
    //sig=1;
    for(int i=1; i<(int)list.size(); i++){
      if(list[i]==list[i-1])return 0;
      if(list[i]<list[i-1]){
        int tmp=list[i-1];
        list[i-1]=list[i];
        list[i]=tmp;
        sig*=-1;
        i--;
        if(i>0)i--;
        continue;
      }
    }
    return sig;
  }
  

  void copy(const Fermion_Index &other){
    list=other.list;
    sig=other.sig;
  }
  Fermion_Index(const Fermion_Index &other){
    copy(other);
  }
  Fermion_Index(const std::vector<int> &list_) : list(list_), sig(1){
  }
  Fermion_Index(): sig(1){
  }
  Fermion_Index(int rank): sig(1){
    list.resize(rank);
    for(int i=0; i<rank; i++)list[i]=i;
  }

};

std::ostream &operator<<(std::ostream &os, const Fermion_Index &fi){
  if(fi.size()<1)return os;
  os<<fi[0];
  for(size_t i=1; i<fi.size(); i++){
    os<<" "<<fi[i];
  }
  return os;
}


Fermion_Index CdaggerC(int i1, int i2, const Fermion_Index &fi_in){
  Fermion_Index fi(fi_in);
  for(int i=0; i<(int)fi_in.size(); i++){
    if(fi_in[i]==i2){
      fi[i]=i1;
      fi.sort();
      return fi;
    }
  }
  fi.sig=0;
  return fi;
}


class Fermion_Enumerator{
public:
  int rank; int nr_states;
  N_Choose_K n_choose_k;
  
  Fermion_Enumerator &set_rank(int r){rank=r; return *this;}
  Fermion_Enumerator &set_nr_states(int s){nr_states=s; return *this;}

  int enumerate_fixed_N(const Fermion_Index &fi, int start_index=0){
  /*
    The 0 index will refer to (0,1,2,..,N-1).
  */ 
    if(fi.size()==0)return 0;
    if(fi.size()==1)return fi[0];
    int subtr=0;
    if(start_index>0){
      subtr=fi[start_index-1]+1;
    } 
    if(fi[start_index]<subtr){
      std::cerr<<"enumerate_fixed_N: unsorted Fermion_Index "<<fi<<"!"<<std::endl;
      exit(1);
    }
    if(start_index==fi.size()-1){
      return fi[start_index]-subtr;
    }

    int res=0;
    for(int i=subtr; i<fi[start_index]; i++){
//      res+=n_choose_k(nr_states-1-subtr, i-subtr+1);
      res+=n_choose_k(nr_states-1-i, fi.size()-1-start_index);
    }
    return res+enumerate_fixed_N(fi, start_index+1);
  }
  int operator()(const Fermion_Index &fi){
    return enumerate_fixed_N(fi);
  }
  int operator()(const Fermion_Index &fi, int nr_states_){
    nr_states=nr_states_;
    return enumerate_fixed_N(fi);
  }

  Fermion_Index get_Fermion_Index(int index){
    Fermion_Index fi(rank);
    if(rank==0)return fi;
    if(rank==1){fi[0]=index;return fi;}

    for(int s=0; s<rank; s++){
      int subtr=0;
      if(s>0){
        subtr=fi[s-1]+1;
      }
      if(s==rank-1){
        fi[s]=index+subtr;
        return fi;
      }
      bool done=false;
      for(int i=subtr; !done; i++){
//        std::cout<<index<<" "<<s<<" "<<subtr<<"->"<<nr_states-(i+1)<<" "<<rank-1-s<<std::endl;
        int nck=n_choose_k(nr_states-1-i, rank-1-s);
//        std::cout<<"nck: "<<nck<<std::endl;
        if(index<nck){
          fi[s]=i;
          done=true;
        }else{
          index-=nck;
        }
      }
    }
    return fi;
  }
  Fermion_Index get_Fermion_Index(int index, int rnk, int nst){
    rank=rnk;
    nr_states=nst;
    return get_Fermion_Index(index);
  }
  int max_index(){
    return n_choose_k(nr_states, rank);
  }
  int max_index(int rnk, int nst){
    rank=rnk; 
    nr_states=nst;
    return n_choose_k(nr_states, rank);
  }

  Fermion_Enumerator(int rank_=0, int nr_states_=0)
   : rank(rank_), nr_states(nr_states_), n_choose_k(nr_states_){
  }

};


#endif
