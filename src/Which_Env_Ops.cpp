#include "PCH.hpp"
#include "Which_Env_Ops.hpp"
#include "Parameters.hpp"
#include "Reader.hpp"
#include "InitialState.hpp"
#include "DummyException.hpp"

namespace ACE{

  void Which_Env_Ops_List::setup(Parameters &param, int sysdim){
    if(sysdim<1){
      Eigen::MatrixXcd initial_rho=InitialState(param);
      sysdim=initial_rho.rows();
    }
    int rows=param.get_nr_rows("add_Env_Op");
    for(int r=0; r<rows; r++){
      std::vector<std::string> row=param.get_row("add_Env_Op", r);
      if(row.size()<2){
        std::cerr<<"Usage: add_Env_Op NR_IF NR_ENV_OP [system_operator]"<<std::endl;
        throw DummyException();
      }
      Which_Env_Ops w;
      w.i=readSizeT(row[0],"add_Env_Op NR_IF");
      w.o=readSizeT(row[1],"add_Env_Op NR_ENV_OP");
      if(row.size()>2){
        w.A=ReadExpression(row[2]);
      }else{
        w.A=Eigen::MatrixXcd::Identity(sysdim, sysdim);
      }
      list.push_back(w);
    }
  }

}//namespace
