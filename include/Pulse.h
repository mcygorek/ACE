#ifndef PULSE_DEFINED_H
#define PULSE_DEFINED_H

#include "Function.h"
#include "Constants.h"
#include "Parameters.h"


/// Pulse envelope: constant
class Pulse_Constant: public ComplexFunction{
public:
  std::complex<double> c;
  double detuning;

  virtual std::complex<double> f(double x)const{ 
    return c*exp(std::complex<double>(0, -detuning*x/Constants::hbar_in_meV_ps)); 
  }
  Pulse_Constant(std::complex<double> c_=0., int detuning_=0) 
   : c(c_), detuning(detuning_){
  }
};


/// Standard Gaussian pulse: area in units of PI
class Pulse_Gauss: public ComplexFunction{
public:
  //area in units of pi; dw in meV; c and FWHM in ps
  double c;
  double FWHM;
  double area;
  double dw;

  virtual std::complex<double> f(double x)const{
    double s=FWHM/(2.*sqrt(2.*log(2.)));
    double y=(x-c)/s;
    if(y*y>1e2)return 0.;
    return exp(-0.5*y*y)/(sqrt(2.*M_PI)*s)*area*M_PI*Constants::hbar_in_meV_ps *exp(std::complex<double>(0,-dw*x/Constants::hbar_in_meV_ps));
  }

  void setup(const std::vector<std::string> &sv){
    if(sv.size()<3){
      std::cerr<<"'pulse_gauss' needs 3 arguments: center, FWHM, area, [detuning]!"<<std::endl;
      exit(1);
    }
    if(sv.size()>3)dw=Reader::readDouble(sv[3],"pulse_gauss: detuning");
///Constants::hbar_in_meV_ps;
    
    c=Reader::readDouble(sv[0],"pulse_gauss: center");
    FWHM=Reader::readDouble(sv[1],"pulse_gauss: FWHM");
    area=Reader::readDouble(sv[2],"pulse_gauss: area");
  }
  Pulse_Gauss(const std::vector<std::string> &sv){
    setup(sv);
  }
  Pulse_Gauss(double c_, double FWHM_, double area_, double dw_=0)
   : c(c_), FWHM(FWHM_), area(area_), dw(dw_){
  }
};


/// Rectangular pulse (time domain): area in units of PI
class Pulse_Rect: public ComplexFunction{
public:
  //area in units of pi; dw in meV; c and FWHM in ps
  double c;
  double width;
  double area;
  double dw;

  virtual std::complex<double> f(double x)const{
    if( fabs(x-c) >  width/2.) return 0.;
    
    return area*M_PI/width*Constants::hbar_in_meV_ps 
        *exp(std::complex<double>(0,-dw*x/Constants::hbar_in_meV_ps));
  }

  void setup(const std::vector<std::string> &sv){
    if(sv.size()<3){
      std::cerr<<"'pulse_rect' needs 3 arguments: center, width, area, [detuning]!"<<std::endl;
      exit(1);
    }
    if(sv.size()>3)dw=Reader::readDouble(sv[3],"pulse_rect: detuning")/Constants::hbar_in_meV_ps;
    
    c=Reader::readDouble(sv[0],"pulse_rect: center");
    width=Reader::readDouble(sv[1],"pulse_rect: width");
    area=Reader::readDouble(sv[2],"pulse_rect: area");
  }
  Pulse_Rect(const std::vector<std::string> &sv){
    setup(sv);
  }
  Pulse_Rect(double c_, double width_, double area_, double dw_=0)
   : c(c_), width(width_), area(area_), dw(dw_){
  }
};




///Rectangle in frequency space -> sinc
class Pulse_rect_freq: public ComplexFunction{
public:
  double c;
  double area;
  double dw;
  double detuning;

  //note: int sin(x)/x = 1/PI  -> normalize
  virtual std::complex<double> f(double x)const{
    double t=(x-c)/Constants::hbar_in_meV_ps;  //in inverse meV;
    double sc=1;
    if(t*t>1e-10){
      sc=sin(dw/2.*t)/(dw/2.*t);
    }
    return dw/2.*area*exp(std::complex<double>(0, -detuning*t))*sc;
  }

  Pulse_rect_freq(double c_, double area_, double dw_, double detuning_)
   : c(c_), area(area_), dw(dw_), detuning(detuning_) {
  }
};



void Pulses_print(const std::vector<ComplexFunctionPtr> & pulses, 
                  const std::string &fname, 
                  double ta, double dt, double te){
  if(fname==""){
    std::cerr<<"Pulses_print: no file specified!"<<std::endl; 
    exit(1); 
  }
  std::ofstream ofs(fname.c_str());

  int Ndiscr=(te-ta+0.5*dt)/dt;
  for(int i=0; i<=Ndiscr; i++){
    double t=ta+i*dt;
    ofs<<t;
    for(size_t j=0; j<pulses.size(); j++){
      std::complex<double> c=pulses[j]->f(t);
      ofs<<" "<<c.real()<<" "<<c.imag();
    }
    ofs<<std::endl;
  }
}

void Pulses_print_sum(const std::vector<ComplexFunctionPtr> & pulses, 
                  const std::string &fname, 
                  double ta, double dt, double te){
  if(fname==""){
    std::cerr<<"Pulses_print: no file specified!"<<std::endl; 
    exit(1); 
  }
  std::ofstream ofs(fname.c_str());

  int Ndiscr=(te-ta+0.5*dt)/dt;
  for(int i=0; i<=Ndiscr; i++){
    double t=ta+i*dt;
    std::complex<double> sum=0;
    for(size_t j=0; j<pulses.size(); j++){
      sum+=pulses[j]->f(t);
    }
    ofs<<t<<" "<<sum.real()<<" "<<sum.imag()<<std::endl;
  }
}


/*
std::vector<ComplexFunctionPtr> Pulses_from_Parameters(Parameters &param){
  std::vector<ComplexFunctionPtr> pulses;

  if(param.is_specified("pulse_gauss")){
    const Parameters_Entry p=param.get("pulse_gauss");
    for(size_t r=0; r<p.size(); r++){
      if(p[r].size()<3){
        std::cerr<<"'pulse_gauss' needs 3 arguments: center, FWHM, area!"<<std::endl;
        exit(1);
      }
      double dw=0;
      if(p[r].size()>3)dw=Reader::readDouble(p[r][3],"pulse_gauss: dw");

      ComplexFunctionPtr pulse=new Pulse_Gauss(
         Reader::readDouble(p[r][0],"pulse_gauss: center"),
         Reader::readDouble(p[r][1],"pulse_gauss: FWHM"),
         Reader::readDouble(p[r][2],"pulse_gauss: area"), dw  );
      
      pulses.push_back(pulse);
    }
  }


  if(param.is_specified("pulse_dichromatic")){
    const Parameters_Entry p=param.get("pulse_dichromatic");
    for(size_t r=0; r<p.size(); r++){
      if(p[r].size()<4){
        std::cerr<<"'pulse_dichromatic' needs 4 [5] arguments: center, FWHM, area, dw, [phi]!"<<std::endl;
        exit(1);
      }
      double phi=0;
      if(p[r].size()>4)phi=Reader::readDouble(p[r][4],"pulse_dichromatic: phi")*M_PI;
      ComplexFunctionPtr pulse_gauss=new Pulse_Gauss(
         Reader::readDouble(p[r][0],"pulse_dichromatic: center"),
         Reader::readDouble(p[r][1],"pulse_dichromatic: FWHM"),
         Reader::readDouble(p[r][2],"pulse_dichromatic: area")  );
      ComplexFunctionPtr pulse=new Pulse_Dichromatic( pulse_gauss, 
         Reader::readDouble(p[r][3],"pulse_dichromatic: dw"), phi );
    
      pulses.push_back(pulse);
    }
  }


//print if advised by parameters:
  if(param.is_specified("print_pulse")){
    double ta=param.get_as_double("ta");
    double dt=param.get_as_double("dt");
    int Nintermediate=param.get_as_double("Nintermediate");
    double te=param.get_as_double("te");
    Pulses_print(pulses, param.get_as_string("print_pulse"),
                 ta, dt/(Nintermediate+1), te);
  }
  if(param.is_specified("print_pulse_sum")){
    double ta=param.get_as_double("ta");
    double dt=param.get_as_double("dt");
    int Nintermediate=param.get_as_double("Nintermediate");
    double te=param.get_as_double("te");
    Pulses_print_sum(pulses, param.get_as_string("print_pulse_sum"),
                 ta, dt/(Nintermediate+1), te);
  }


  return pulses;
}
*/


#include "Pulse_Dichromatic.h"
#endif
