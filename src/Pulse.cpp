#include "Pulse.hpp"
#include "ComplexFunction_Interpolate.hpp"
#include "Function.hpp"
#include "Constants.hpp"
#include "Parameters.hpp"
#include "Reader.hpp"

namespace ACE{

/// Pulse envelope: constant
  std::complex<double> Pulse_Constant::f(double x)const{ 
    return c*exp(std::complex<double>(0, -detuning*x/hbar_in_meV_ps)); 
  }


/// Standard Gaussian pulse: area in units of PI
  std::complex<double> Pulse_Gauss::f(double x)const{
    double s=FWHM/(2.*sqrt(2.*log(2.)));
    double y=(x-c)/s;
    if(y*y>1e2)return 0.;
    return area * exp(-0.5*y*y)/(sqrt(2.*M_PI)*s)
          * exp(std::complex<double>(0,-dE*x/hbar_in_meV_ps));
  }

  void Pulse_Gauss::setup(const std::vector<std::string> &sv){
    if(sv.size()<3){
      std::cerr<<"'pulse_gauss' needs 3 arguments: center, FWHM, area, [detuning]!"<<std::endl;
      exit(1);
    }
    dE=0.;
    if(sv.size()>3)dE=readDouble(sv[3],"pulse_gauss: detuning");
    
    c=readDouble(sv[0],"pulse_gauss: center");
    FWHM=readDouble(sv[1],"pulse_gauss: FWHM");
    area=readDouble(sv[2],"pulse_gauss: area");
  }

  std::complex<double> Pulse_SmoothRect::f(double x)const{
    double a=0.;
    if(t_rise<=0.){
      a=(x<t_on || x>t_off) ? 0. : 1.; 
    }else{
      if(-(x-t_on) > 25.*t_rise || -(t_off-x) > 25.*t_rise) return 0.; 
      a=1./( (1.+exp(-(x-t_on)/t_rise))*(1.+exp(-(t_off-x)/t_rise)));
    }
    
    return scale*a*exp(std::complex<double>(0,-detuning*x/hbar_in_meV_ps));
  }


/// Rectangular pulse (time domain): area in units of PI
  std::complex<double> Pulse_Rect::f(double x)const{
    if( fabs(x-c) >  width/2.) return 0.;
    
    return (area/width)
        *exp(std::complex<double>(0,-dE*x/hbar_in_meV_ps));
  }

  void Pulse_Rect::setup(const std::vector<std::string> &sv){
    if(sv.size()<3){
      std::cerr<<"'pulse_rect' needs 3 arguments: center, width, area, [detuning]!"<<std::endl;
      exit(1);
    }
    if(sv.size()>3)dE=readDouble(sv[3],"pulse_rect: detuning");
    
    c=readDouble(sv[0],"pulse_rect: center");
    width=readDouble(sv[1],"pulse_rect: width");
    area=readDouble(sv[2],"pulse_rect: area");
  }

///Rectangle in frequency space -> sinc

  //note: int sin(x)/x = 1/PI  -> normalize
  std::complex<double> Pulse_rect_freq::f(double x)const{
    double t=(x-c);  //in inverse meV;
    double dw=dE/hbar_in_meV_ps;
    double sc=1;
    if(t*t>1e-10){
      sc=sin(dw/2.*t)/(dw/2.*t);
    }
    return dw/2.*area*exp(std::complex<double>(0, -detuning*t/hbar_in_meV_ps))*sc;
  }

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
      if(p[r].size()>3)dw=readDouble(p[r][3],"pulse_gauss: dw");

      ComplexFunctionPtr pulse=new Pulse_Gauss(
         readDouble(p[r][0],"pulse_gauss: center"),
         readDouble(p[r][1],"pulse_gauss: FWHM"),
         readDouble(p[r][2],"pulse_gauss: area"), dw  );
      
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
      if(p[r].size()>4)phi=readDouble(p[r][4],"pulse_dichromatic: phi")*M_PI;
      ComplexFunctionPtr pulse_gauss=new Pulse_Gauss(
         readDouble(p[r][0],"pulse_dichromatic: center"),
         readDouble(p[r][1],"pulse_dichromatic: FWHM"),
         readDouble(p[r][2],"pulse_dichromatic: area")  );
      ComplexFunctionPtr pulse=new Pulse_Dichromatic( pulse_gauss, 
         readDouble(p[r][3],"pulse_dichromatic: dw"), phi );
    
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


  std::complex<double> Pulse_Dichromatic_Rect::f(double x)const{
    std::complex<double> res=0.;
    for(size_t i=0; i<pulses.size(); i++){
      res+=pulses[i]->f(x);
    }
    return res;
  }

  Pulse_Dichromatic_Rect::Pulse_Dichromatic_Rect(
      double area, double C, double dw, double DW, 
      const std::string &pulse_type, double detuning, double center){ 

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
      double convert=hbar_in_meV_ps*8.*log(2.);

//std::cout<<"eff_area1: "<<eff_area1<<" "<<detuning-DW/2.<<std::endl;
//std::cout<<"eff_area2: "<<eff_area2<<" "<<detuning+DW/2.<<std::endl;
      if(C<1.-1e-8)pulses.push_back(std::make_shared<Pulse_Gauss>(center, convert/dw1, eff_area1, detuning-DW/2.));
      if(C>-1.+1e-8)pulses.push_back(std::make_shared<Pulse_Gauss>(center, convert/dw2, eff_area2, detuning+DW/2.));


    }else if(pulse_type=="gauss_height"){
      double convert=hbar_in_meV_ps*8.*log(2.);
      if(C<0){
        eff_area1=area*sqrt(2*1./(1.+f*f));
        eff_area2=area*sqrt(2*f*f/(1.+f*f));
      }else{
        eff_area2=area*sqrt(2*1./(1.+f*f));
        eff_area1=area*sqrt(2*f*f/(1.+f*f));
      }
      
      if(C<1.-1e-8)pulses.push_back(std::make_shared<Pulse_Gauss>(center, convert/dw, eff_area1, detuning-DW/2.));
      if(C>-1.+1e-8)pulses.push_back(std::make_shared<Pulse_Gauss>(center, convert/dw, eff_area2, detuning+DW/2.));


    }else if(pulse_type=="asym"){
      double convert=hbar_in_meV_ps*8.*log(2.);

      if(C<1.-1e-8){
        pulses.push_back(std::make_shared<Pulse_Gauss>(center, convert/dw1, eff_area1, detuning-DW/2.));
        pulses.push_back(std::make_shared<Pulse_Gauss>(center, convert/(dw1/4.), 0.5*eff_area1, detuning-DW/2.+dw1/2.));
      }
      if(C>-1.+1e-8){
        pulses.push_back(std::make_shared<Pulse_Gauss>(center, convert/dw2, eff_area2, detuning+DW/2.));
        pulses.push_back(std::make_shared<Pulse_Gauss>(center, convert/(dw2/4.), 0.5*eff_area2, detuning+DW/2.-dw2/2.));
      }


    }else if(pulse_type=="gauss_cut"){
      double convert=hbar_in_meV_ps*8.*log(2.);

      pulses.push_back(std::make_shared<Pulse_Gauss>(center, convert/DW, area, detuning));
      pulses.push_back(std::make_shared<Pulse_rect_freq>(center, -area, dw, detuning));


    }else{
      if(C<1.-1e-8)pulses.push_back(std::make_shared<Pulse_rect_freq>(center, eff_area1, dw1, detuning-DW/2.));
      if(C>-1.+1e-8)pulses.push_back(std::make_shared<Pulse_rect_freq>(center, eff_area2, dw2, detuning+DW/2.));
    }
  }



  ComplexFunctionPtr Pulse_from_data(const std::pair<std::vector<double>,std::vector<std::complex<double> > > & shape){
    ComplexFunctionPtr ptr=std::make_shared<ComplexFunction_Interpolate>(); 
    ComplexFunction_Interpolate* p = static_cast<ComplexFunction_Interpolate*>(ptr.get()); 
    int length=shape.first.size(); 
    p->val.resize(length); 
    for(int i=0; i<length; i++){ 
      p->val[i].first=shape.first[i]; 
      p->val[i].second= (i<shape.second.size() ? shape.second[i] : 0.); 
    }
    return ptr;
  }
}//namespace
