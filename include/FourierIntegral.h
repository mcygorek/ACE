#ifndef FOURIER_INTEGRAL_DEFINED_H
#define FOURIER_INTEGRAL_DEFINED_H

#include <fstream>
#include "continuousFFT.h"
#include "Simulation_Results.h"



class FourierIntegral{
public:
  bool was_calculated;
  int Ntruncate;

  Simulation_Results results;
  double wa, dw;
  
 
  void initialize(){
    was_calculated=false;
    Ntruncate=-1;
  }
  void calculate(const Simulation_Results &in, int col, bool forward=true){
    std::vector<std::complex<double> > vout;
    if(in.list.size()<2){
      std::cerr<<"FourierIntegral::calculate: in.list.size()<2 !"<<std::endl;
      exit(1);
    }

    double ta=in.list[0].first;
    double dt=in.list[1].first-in.list[0].first;
    int NN=in.list.size()*10;
    if(Ntruncate>0)NN=Ntruncate;
    if(NN<in.list.size()+1)NN=in.list.size()+1;
    wa=continuousFFT_center_wa(dt, NN);
    dw=continuousFFT_dw(dt, NN);

    std::vector<std::complex<double> > v(in.list.size());
    for(size_t i=0; i<in.list.size(); i++){
      if(in.list[i].second.size()<=col){
        std::cerr<<"continuousFFT vin["<<i<<"].second.size() <= col!"<<std::endl;
        exit(1);
      }
      v[i]=in.list[i].second[col];
    }
    continuousFFT(v, vout, ta, dt, NN, wa, &dw, forward);

//    continuousFFT(in.list, col, vout, ta, dt, NN, wa);

    results.clear();
    results.resize(vout.size());
    for(size_t i=0; i<vout.size(); i++){
      double w=wa+i*dw;
      std::vector<std::complex<double> > v(1., vout[i]);
      results[i]=std::make_pair(w, v);
    } 


    was_calculated=true;
  }
 
  void print(const std::string &fname){
    if(!was_calculated){
      std::cerr<<"FourierIntegral: ERROR: printing without having calculated results!"<<std::endl;
      exit(1);
    }
    results.print(fname);
  }
  
  FourierIntegral(){
    initialize();
  }
  FourierIntegral(const Simulation_Results &in, int col){ 
    initialize();
    calculate(in, col);
  }
};

#endif 
