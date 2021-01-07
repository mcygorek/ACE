#ifndef PULSE_PRINTER_DEFINED_H
#define PULSE_PRINTER_DEFINED_H

#include "Parameters.h"
#include "FreePropagator.h"
#include "slowFT.h"


void Pulse_Printer(Parameters &param, FreePropagator &fprop){
  double ta=param.get_as_double("ta",0.);
  double te=param.get_as_double("te",1.);
  double dt=param.get_as_double("dt",0.1);
  int Nintermediate=param.get_as_size_t("Nintermediate",0);

  if(param.is_specified("print_pulse_FT")){
    std::string errmsg="Usage: print_pulse_FT FILENAME E_low E_high";// Ndiscr";
    std::vector<std::string> sv=param.get_row("print_pulse_FT",0);
    if(sv.size()<1){std::cerr<<errmsg<<std::endl; exit(1);}

    double ftlim_a, ftlim_b;
    if(sv.size()>1){ftlim_a=Reader::readDouble(sv[1], "print_pulse_FT: E_low");}
    else{ ftlim_a=-3.; }
    if(sv.size()>2){ftlim_b=Reader::readDouble(sv[2], "print_pulse_FT: E_high");}
    else{ ftlim_b=3.; }

    int Nsteps=((te-ta)/dt)*(Nintermediate+1)+1;
    double dt2=dt/(Nintermediate+1);
    std::vector<std::complex<double> > pulses_interpol(Nsteps,0.);
    for(int i=0; i<Nsteps; i++){
      double t=ta+i*dt2;
      for(size_t j=0; j<fprop.timedep_H.size(); j++){
        pulses_interpol[i]+=fprop.timedep_H[j].first->f(t);
      }
    }
    print_slowFT_meV(sv[0], pulses_interpol, ftlim_a, ftlim_b, 1000, ta, dt2, 1, 100);
  }
}



#endif
