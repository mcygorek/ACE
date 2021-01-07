#ifndef INITIAL_STATE_DEFINED_H
#define INITIAL_STATE_DEFINED_H

#include "Parameters.h"
#include "Operators.h"
#include "ReadOperator.h"

class InitialState{
public:

  Eigen::MatrixXcd rho;

  operator Eigen::MatrixXcd&(){
    return rho;
  }
  void initialize(){
    rho=Operators2x2::ketbra(0,0);
  }
  void setup(Parameters &param){
    initialize();
    int MTLS=param.get_as_size_t("MTLS",0);
    std::string str=param.get_as_string("initial","");
    if(str==""||str=="TLS_GS"){
      rho=Operators2x2::ketbra(0,0);
      std::cout<<"Using initial state: 2LS ground state"<<std::endl;
    }else if(str=="TLS_EX"||str=="EX"){
      rho=Operators2x2::ketbra(1,1);
      std::cout<<"Using initial state: 2LS excited"<<std::endl;
    }else if(str=="TLS_COH"){
      rho=Eigen::MatrixXcd(2,2);
      rho<< 0.5, 0.5, 0.5, 0.5;
      std::cout<<"Using initial state: 2LS coherent"<<std::endl;
    }else if(str=="MTLS_EX_ALL"){
      rho=Eigen::MatrixXcd::Zero(loop_pow(2, MTLS), loop_pow(2, MTLS));
      rho(rho.rows()-1,rho.cols()-1)=1.;
      std::cout<<"Using initial state: "<<MTLS<<" 2LSs excited state"<<std::endl;
    }else if(str=="MTLS_GS_ALL"||MTLS>0){
      rho=Eigen::MatrixXcd::Zero(loop_pow(2, MTLS), loop_pow(2, MTLS));
      rho(0,0)=1.;
      std::cout<<"Using initial state: "<<MTLS<<" 2LSs ground state"<<std::endl;
    }else{
//      std::string singlestr=param.get_as_single_string("initial");
//      std::cout<<"Using explicit initial state: '"<<singlestr<<"'"<<std::endl;
////      rho=ReadOperator(singlestr);
//      rho=ReadExpression(singlestr);
      std::cout<<"Using explicit initial state: '"<<str<<"'"<<std::endl;
      rho=param.get_as_operator("initial");
      std::cout<<rho<<std::endl;
    }
  }
  InitialState(Parameters &param){
    setup(param);
  }
  InitialState(){
    initialize();
  }
};

#endif
