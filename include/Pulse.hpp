#ifndef PULSE_DEFINED_H
#define PULSE_DEFINED_H

#include "Function.hpp"
#include <vector>
#include <iosfwd>

namespace ACE{
class Parameters;

/// Pulse envelope: constant
class Pulse_Constant: public ComplexFunction{
public:
  std::complex<double> c;
  double detuning;

  virtual std::complex<double> f(double x)const;

  inline Pulse_Constant(std::complex<double> c_=0., int detuning_=0) 
   : c(c_), detuning(detuning_){
  }
};


/// Standard Gaussian pulse: area in units of PI
class Pulse_Gauss: public ComplexFunction{
public:
  //area in units of 1 (as opposed to pi); dE in meV; c and FWHM in ps
  double c;
  double FWHM;
  double area;
  double dE;

  virtual std::complex<double> f(double x)const;

  void setup(const std::vector<std::string> &sv);

  inline Pulse_Gauss(const std::vector<std::string> &sv): c(0.),FWHM(1.){
    setup(sv);
  }
  inline Pulse_Gauss(double c_, double FWHM_, double area_, double dE_=0)
   : c(c_), FWHM(FWHM_), area(area_), dE(dE_){
  }
};


class Pulse_SmoothRect: public ComplexFunction{
public:
  // detuning in meV; t_on, t_off, t_rise in ps
  double t_on;
  double t_off;
  double t_rise;
  double scale;
  double detuning;

  virtual std::complex<double> f(double x)const;

  Pulse_SmoothRect(double t_on_, double t_off_, double t_rise_=0., 
                   double scale_=1., double detuning_=0.)
   : t_on(t_on_), t_off(t_off_), t_rise(t_rise_), scale(scale_), detuning(detuning_){
  }
};


/// Rectangular pulse (time domain): area in units of PI
class Pulse_Rect: public ComplexFunction{
public:
  //area ; dE in meV; c and FWHM in ps
  double c;
  double width;
  double area;
  double dE;

  virtual std::complex<double> f(double x)const;

  void setup(const std::vector<std::string> &sv);
  
  inline Pulse_Rect(const std::vector<std::string> &sv){
    setup(sv);
  }
  inline Pulse_Rect(double c_, double width_, double area_, double dE_=0)
   : c(c_), width(width_), area(area_), dE(dE_){
  }
};


///Rectangle in frequency space -> sinc
class Pulse_rect_freq: public ComplexFunction{
public:
  double c;
  double area;
  double dE;
  double detuning;

  //note: int sin(x)/x = 1/PI  -> normalize
  virtual std::complex<double> f(double x)const;

  inline Pulse_rect_freq(double c_, double area_, double dE_, double detuning_)
   : c(c_), area(area_), dE(dE_), detuning(detuning_) {
  }
};

class Pulse_Dichromatic_Rect: public ComplexFunction{
public:
  std::vector<ComplexFunctionPtr> pulses;

  virtual std::complex<double> f(double x)const;

  Pulse_Dichromatic_Rect(double area, double C, double dw, double DW, const std::string &pulse_type="rect", double detuning=0, double center=0.);
  
};



void Pulses_print(const std::vector<ComplexFunctionPtr> & pulses, 
                  const std::string &fname, 
                  double ta, double dt, double te);

void Pulses_print_sum(const std::vector<ComplexFunctionPtr> & pulses, 
                  const std::string &fname, 
                  double ta, double dt, double te);


}//namespace
#endif
