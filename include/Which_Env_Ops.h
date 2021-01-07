#ifndef PRINT_WHICH_ENV_OPS_DEFINED_H
#define PRINT_WHICH_ENV_OPS_DEFINED_H

#include "Parameters.h"



struct Which_Env_Ops{
  int i; //which IF
  int o; //which env_op
  Eigen::MatrixXcd A; //which system operator
};

class Which_Env_Ops_List{
public:
  std::vector<Which_Env_Ops> list;

  Which_Env_Ops & operator[] (size_t i){return list[i];}
  size_t size()const{return list.size();}
  
  void setup(Parameters &param){
    int rows=param.get_nr_rows("add_Env_Op");
    for(int r=0; r<rows; r++){
      std::vector<std::string> row=param.get_row("add_Env_Op", r);
      if(row.size()<2){
        std::cerr<<"Usage: add_Env_Op NR_IF NR_ENV_OP [system_operator]"<<std::endl;
        exit(1);
      }
      Which_Env_Ops w;
      w.i=Reader::readSizeT(row[0],"add_Env_Op NR_IF");
      w.o=Reader::readSizeT(row[1],"add_Env_Op NR_ENV_OP");
      if(row.size()>2){
        w.A=ReadExpression(row[2]);
      }
      list.push_back(w);
    }
  }

};

#endif
