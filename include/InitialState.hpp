#pragma once
#ifndef INITIAL_STATE_DEFINED_H
#define INITIAL_STATE_DEFINED_H

#include "Parameters.hpp"
#include "Operators.hpp"
#include "ReadTemperature.hpp"
#include "LiouvilleTools.hpp"

namespace ACE{

class InitialState{
public:

  Eigen::MatrixXcd rho;

  inline operator Eigen::MatrixXcd&(){
    return rho;
  }
  inline void initialize(){
    rho=Operators(2).ketbra(0,0);
  }
  void setup(Parameters &param);
  
  inline InitialState(Parameters &param){
    setup(param);
  }
  inline InitialState(){
    initialize();
  }
};

}//namespace
#endif
