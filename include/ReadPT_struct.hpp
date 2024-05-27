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
