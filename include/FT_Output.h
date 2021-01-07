#ifndef FT_OUTPUT_DEFINED_H
#define FT_OUTPUT_DEFINED_H

#include "FT_Parameters.h"
#include "Parameters.h"
#include "Simulation_Results.h"


class FT_Output{
public:
  std::vector<int> FT_columns;
  std::pair<bool, double> FT_ta;
  std::string FT_file;
  FT_Parameters FT_param;

  void setup(Parameters &param){
    if(param.is_specified("FT_columns")){
      std::vector<std::string> FT_str=param.get_all_strings("FT_columns");
      for(size_t i=0; i<FT_str.size(); i++){
        FT_columns.push_back(Reader::readSizeT(FT_str[i], "FT_columns"));
      }
    }
    if(param.is_specified("FT_ta")){
      FT_ta.first=true;
      FT_ta.second=param.get_as_double("FT_ta");
    }else{
      FT_ta=std::make_pair(false, 0.);
    }
    FT_file="";
    if(param.is_specified("FT_parameters")){
      std::vector<std::string> row=param.get_row("FT_parameters",0);
      if(row.size()<4){
        std::cerr<<"'FT_parameters': FILENAME Emin Emax Ndiscr [Nsubdiv]"<<std::endl;
        exit(1);
      }
      FT_file=row[0];
      FT_param.wa=Reader::readDouble(row[1],"FT_parameters: wa")/Constants::hbar_in_meV_ps;
      FT_param.we=Reader::readDouble(row[2],"FT_parameters: we")/Constants::hbar_in_meV_ps;
      FT_param.Ndiscr=Reader::readSizeT(row[3],"FT_parameters: Ndiscr");
      if(row.size()>4){
        FT_param.Nsubdiv=Reader::readSizeT(row[4],"FT_parameters: Nsubdiv");
      }
    }
  }

  void print(const Simulation_Results &res)const{
    if(FT_file=="")return; //do nothing
    if(FT_columns.size()<1)return; //do nothing

    if(res.size()<2){
      std::cerr<<"Error: FT_Output: Need at least two data points!"<<std::endl;
      exit(1);
    }
    for(size_t i=0; i<FT_columns.size(); i++){
      if(FT_columns[i]>=res[0].second.size()){
        std::cerr<<"Error: FT_Output: column index "<<FT_columns[i]<<" too large!"<<std::endl;
        exit(1);
      }
    }
 
    FT_Parameters FT_param(this->FT_param);

    FT_param.ta=res[0].first;
    FT_param.dt=res[1].first-FT_param.ta;
    int ignore_first=0;   
 
    if(FT_ta.first){
      ignore_first=(FT_ta.second-FT_param.ta)/FT_param.dt;
      FT_param.ta=0.;//FT_ta.second;
      if(((int)res.size())-ignore_first<2){
        std::cerr<<"Error: FT_Output: !((int)FT_columns.size())-ignore_first<2"<<std::endl;
        std::cerr<<"ignore_first: "<<ignore_first<<std::endl;
        std::cerr<<"FT_param.ta: "<<FT_param.ta<<std::endl;
        std::cerr<<"FT_ta.second: "<<FT_ta.second<<std::endl;
        std::cerr<<"FT_param.dt: "<<FT_param.dt<<std::endl;
        exit(1);
      }
    }
    
    
    Simulation_Results FT_res;
    FT_res.resize(FT_param.Ndiscr+1);

    //setup energy discretization
    for(size_t i=0; i<FT_res.size(); i++){
      double lambda=i/((double)FT_param.Ndiscr);
      FT_res[i].first=Constants::hbar_in_meV_ps*(FT_param.wa*(1.-lambda)+FT_param.we*lambda);
    }

    //FT column by column:
    for(size_t c=0; c<FT_columns.size(); c++){
      std::vector<std::complex<double> > in;
      for(int l=ignore_first; l<res.size(); l++){
        in.push_back(res[l].second[FT_columns[c]]);
      }
      std::vector<std::complex<double> > sp=slowFT(in, FT_param);
      for(size_t i=0; i<FT_res.size(); i++){
        FT_res[i].second.push_back(sp[i]);
      }
    }

    FT_res.print(FT_file);
  }

  FT_Output(Parameters &param){
    setup(param);
  }
  FT_Output(){
    Parameters param;
    setup(param);
  }
};


#endif
