#pragma once
#ifndef FT_OUTPUT_DEFINED_H
#define FT_OUTPUT_DEFINED_H

#include "FT_Parameters.hpp"
#include "Parameters.hpp"
#include "Simulation_Results.hpp"

namespace ACE{

class FT_Output{
public:
  std::vector<int> FT_columns;
  std::pair<bool, double> FT_ta;
  std::string FT_file;
  FT_Parameters FT_param;

  void setup(Parameters &param);

  void print(const Simulation_Results &res)const;

  inline FT_Output(Parameters &param){
    setup(param);
  }
  inline FT_Output(){
    Parameters param;
    setup(param);
  }
};

}//namespace
#endif
