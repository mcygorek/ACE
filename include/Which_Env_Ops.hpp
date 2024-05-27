#pragma once
#ifndef PRINT_WHICH_ENV_OPS_DEFINED_H
#define PRINT_WHICH_ENV_OPS_DEFINED_H

#include "Parameters.hpp"

namespace ACE{

struct Which_Env_Ops{
  int i; //which IF
  int o; //which env_op
  Eigen::MatrixXcd A; //which system operator
};

class Which_Env_Ops_List{
public:
  std::vector<Which_Env_Ops> list;

  inline Which_Env_Ops & operator[] (size_t i){return list[i];}
  inline int get_dim()const{if(list.size()>0)return list[0].A.rows(); else return 0;}
  inline size_t size()const{return list.size();}
  inline std::vector<Which_Env_Ops>::const_iterator begin()const{return list.begin();}
  inline std::vector<Which_Env_Ops>::const_iterator end()const{return list.end();}
  
  void setup(Parameters &param, int sysdim=-1);
};

}//namespace
#endif
