#pragma once
#ifndef SIMULATION_RESULTS_DEFINED_H
#define SIMULATION_RESULTS_DEFINED_H
#include "FT_Parameters.hpp"
#include "Output_Ops.hpp"
#include <fstream>

namespace ACE{

typedef std::pair<double, std::vector<std::complex<double> > > Simulation_Results_Entry;


class Simulation_Results{
public: 
  std::vector<Simulation_Results_Entry> list;

  inline size_t size()const{return list.size();}
  inline void clear(){return list.clear();}
  inline void resize(size_t s){return list.resize(s);}
  inline Simulation_Results_Entry & operator[](size_t i){return list[i];}
  inline const Simulation_Results_Entry & operator[](size_t i)const {return list[i];}
  inline void push_back(const Simulation_Results_Entry & r){list.push_back(r);}
  inline Simulation_Results_Entry & back(){return list.back();}
  inline const Simulation_Results_Entry & back()const{return list.back();}

  //add contributions from multiple runs
  Simulation_Results & combine(const Simulation_Results &other);

  void print_single(std::ostream &ofs, int i)const;
  
  void print(std::ostream &ofs)const;
  
  inline void print(const std::string &fname)const{
    std::ofstream ofs(fname.c_str());
    print(ofs);
  }
  void read(const std::string &fname);

  void set(int step, double t, const Output_Ops & output_Op, const Eigen::MatrixXcd &rho, const Eigen::MatrixXcd *Hamil=NULL);

  void add_back(int step, const std::complex<double> &c);
  
  void add_back_nan(int step, int how_many=1);

  inline Simulation_Results(){}
  inline Simulation_Results(size_t N):list(N){}
  inline Simulation_Results(const std::string &str){
    read(str);
  }
};

}//namespace
#endif
