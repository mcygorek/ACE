#ifndef PULSE_DICHROMATIC_DEFINED_H
#define PULSE_DICHROMATIC_DEFINED_H

#include "Pulse.h"


class Pulse_Dichromatic_Rect: public ComplexFunction{
public:
  std::vector<ComplexFunctionPtr> pulses;

  virtual std::complex<double> f(double x)const{
    std::complex<double> res=0.;
    for(size_t i=0; i<pulses.size(); i++){
      res+=pulses[i]->f(x);
    }
    return res;
  }

  Pulse_Dichromatic_Rect(double area, double C, double dw, double DW, const std::string &pulse_type="rect", double detuning=0, double center=0.){ 


    double eff_area1=area;
    double eff_area2=area;
    double dw1=dw;
    double dw2=dw;
    double absC=fabs(C);
    double f=(1.-absC)/(1.+absC);
    if(C<0){
      dw2=f*dw;
    }else{
      dw1=f*dw;
    }
    eff_area1*=sqrt((2.*dw)/(dw1+dw2));
    eff_area2*=sqrt((2.*dw)/(dw1+dw2));

//std::cout<<"C: "<<C<<" f: "<<f<<" dw1: "<<dw1<<" dw2: "<<dw2;
//std::cout<<" eff_area1: "<<eff_area1<<" eff_area2: "<<eff_area2<<std::endl;

    if(pulse_type=="gauss"){
      //note: FWHM=2.*sqrt(2.*log(2.))*sigma;  sigma_FT=1/sigma
      // => FWHM_FT=2.*sqrt(2.*log(2.))*sigma_FT=[2.*sqrt(2.*log(2.))]^2/FWHM
      double convert=Constants::hbar_in_meV_ps*8.*log(2.);

//std::cout<<"eff_area1: "<<eff_area1<<" "<<detuning-DW/2.<<std::endl;
//std::cout<<"eff_area2: "<<eff_area2<<" "<<detuning+DW/2.<<std::endl;
      if(C<1.-1e-8)pulses.push_back(new Pulse_Gauss(center, convert/dw1, eff_area1, detuning-DW/2.));
      if(C>-1.+1e-8)pulses.push_back(new Pulse_Gauss(center, convert/dw2, eff_area2, detuning+DW/2.));


    }else if(pulse_type=="gauss_height"){
      double convert=Constants::hbar_in_meV_ps*8.*log(2.);
      if(C<0){
        eff_area1=area*sqrt(2*1./(1.+f*f));
        eff_area2=area*sqrt(2*f*f/(1.+f*f));
      }else{
        eff_area2=area*sqrt(2*1./(1.+f*f));
        eff_area1=area*sqrt(2*f*f/(1.+f*f));
      }
      
      if(C<1.-1e-8)pulses.push_back(new Pulse_Gauss(center, convert/dw, eff_area1, detuning-DW/2.));
      if(C>-1.+1e-8)pulses.push_back(new Pulse_Gauss(center, convert/dw, eff_area2, detuning+DW/2.));


    }else if(pulse_type=="asym"){
      double convert=Constants::hbar_in_meV_ps*8.*log(2.);

      if(C<1.-1e-8){
        pulses.push_back(new Pulse_Gauss(center, convert/dw1, eff_area1, detuning-DW/2.));
        pulses.push_back(new Pulse_Gauss(center, convert/(dw1/4.), 0.5*eff_area1, detuning-DW/2.+dw1/2.));
      }
      if(C>-1.+1e-8){
        pulses.push_back(new Pulse_Gauss(center, convert/dw2, eff_area2, detuning+DW/2.));
        pulses.push_back(new Pulse_Gauss(center, convert/(dw2/4.), 0.5*eff_area2, detuning+DW/2.-dw2/2.));
      }


    }else if(pulse_type=="gauss_cut"){
      double convert=Constants::hbar_in_meV_ps*8.*log(2.);

      pulses.push_back(new Pulse_Gauss(center, convert/DW, area, detuning));
      pulses.push_back(new Pulse_rect_freq(center, -area, dw, detuning));


    }else{
      if(C<1.-1e-8)pulses.push_back(new Pulse_rect_freq(center, eff_area1, dw1, detuning-DW/2.));
      if(C>-1.+1e-8)pulses.push_back(new Pulse_rect_freq(center, eff_area2, dw2, detuning+DW/2.));
    }

  }

};





#endif
