#ifndef ACE_READ_PT_STRUCT_DEFINED_H
#define ACE_READ_PT_STRUCT_DEFINED_H

#include "Parameters.hpp"

namespace ACE{

struct ReadPT_struct{
  std::string fname; 
  int expand_front; 
  int expand_back;
  inline bool operator!()const{ return fname=="";}
  inline bool have_to_expand()const{ return (expand_front>1 || expand_back>1);}


  inline int get_front_dim()const{
    if(expand_front>1)return expand_front; 
    else return 1;
  }
  inline int get_back_dim()const{
    if(expand_back>1)return expand_back; 
    else return 1;
  }
  inline int get_total_dim(int Nsub)const{  //total Hilbert space size
    int Ntot=Nsub;
    if(expand_front>1)Ntot*=expand_front;
    if(expand_back>1)Ntot*=expand_back;
    return Ntot;
  }
  
  int get_reduced_index_H(int I, int Nsub)const;  //Hilbert space index
  int get_reduced_index_L(int I, int Nsub)const;  //Liouville space index
  int replace_subsystem_index_H(int I, int j, int Nsub)const; 
  int replace_subsystem_index_L(int I, int j, int Nsub)const; 
  
  void print_info(std::ostream &ofs=std::cout)const;
  //typically: par_name="read_PT", expand_name="read_PT_expand"
  void setup(Parameters &param, const std::string &par_name, const std::string & expand_name);

  ReadPT_struct(const std::string fname_="", int front_=1, int back_=1)
    : fname(fname_), expand_front(front_), expand_back(back_){
  }
  ReadPT_struct(Parameters &param, const std::string &par_name, const std::string & expand_name){
    setup(param, par_name, expand_name);
  }
};


}//namespace
#endif
