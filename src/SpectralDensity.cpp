#include "SpectralDensity.hpp"
#include "Function.hpp"
#include <cstdlib>
#include <iostream>
#include "Constants.hpp"
#include "Parameters.hpp"

namespace ACE{

  std::string SpectralDensity_QD::add_name(const std::string &mpgname, const std::string &str)const{
    if(mpgname=="")return str;
    return std::string(mpgname+"_"+str);
  } 
 
  void SpectralDensity_QD::setup(Parameters &param, const std::string &prefix){
    D_e = param.get_as_double(add_name(prefix,"D_e"),7.0) * 1000;
    D_h = param.get_as_double(add_name(prefix,"D_h"),-3.5) * 1000;
    mass_density = param.get_as_double(add_name(prefix,"mass_density"),5370) *
                                             kg_to_meV_ps2_by_nm2 * 1.0e-27;
    c_s = param.get_as_double(add_name(prefix,"c_s"),5110) * 1e-3;
    a_e = param.get_as_double(add_name(prefix,"a_e"),4);
    a_h = param.get_as_double(add_name(prefix,"a_h"),a_e/1.15);
  }

  double SpectralDensity_QD::f(double x)const{
    double x2=x*x;
    double c2=c_s*c_s;
    double arg=D_e*exp(-x2*a_e*a_e/(4.*c2))-D_h*exp(-x2*a_h*a_h/(4.*c2));
    return x2*x/(4.*M_PI*M_PI*mass_density*hbar_in_meV_ps*c2*c2*c_s)*arg*arg;
  }

  SpectralDensity_QD::SpectralDensity_QD(){
    Parameters param;
    setup(param);
  }
  SpectralDensity_QD::SpectralDensity_QD(Parameters &param, const std::string &prefix){
    setup(param, prefix);
  }


/** (sub-/super-)ohmic spectral density  */
  std::string SpectralDensity_sohmic::add_name(const std::string &mpgname, const std::string &str)const{
    if(mpgname=="")return str;
    return std::string(mpgname+"_"+str);
  } 
 
  void SpectralDensity_sohmic::setup(Parameters &param, const std::string &prefix){
    s = param.get_as_double(add_name(prefix,"s"),1);
    alpha = param.get_as_double_check(add_name(prefix,"alpha"));
    omega_c = param.get_as_double(add_name(prefix,"omega_c"));
    div_by_cutoff = param.get_as_bool(add_name(prefix,"div_by_cutoff"),false);

    std::string cstr=param.get_as_string(add_name(prefix,"cutoff"),"none");
    if(cstr=="Exp"||cstr=="exp"){ cutoff_mode=EXP; }
    else if(cstr=="Gauss"||cstr=="gauss"){ cutoff_mode=GAUSS; }
    else if(cstr=="Drude"||cstr=="drude"){ cutoff_mode=DRUDE; }
    else if(cstr=="Hard"||cstr=="hard"){ cutoff_mode=HARD; }
    else if(cstr=="None"||cstr=="none"){ cutoff_mode=NONE; }
    else{
      std::cerr<<"SpectralDensity_sohmic: unknown cutoff mode '"<<cstr<<"'!"<<std::endl;
      exit(1);
    }
    std::cout<<"Using ("<<s<<"-)ohmic spectral density with cut-off '"<<cstr<<"', alpha="<<alpha<<", omega_c="<<omega_c<<std::endl;
  }

  double SpectralDensity_sohmic::f(double x)const{
    double fco=1.;
    if(cutoff_mode==GAUSS){ fco=exp(- (x/omega_c)*(x/omega_c)); }
    else if(cutoff_mode==DRUDE){ fco=1./(1.+(x/omega_c)*(x/omega_c)); }
    else if(cutoff_mode==EXP){ fco=exp(-fabs(x)/omega_c); }
    else if(cutoff_mode==HARD){ fco=(x<=omega_c ? 1 : 0); }
     
    if(fabs(x)<1e-8){
      return 0;
    }else if(div_by_cutoff){
      return alpha*x*pow(x/omega_c,s-1.)*fco;
    }else{
      return alpha*pow(x,s)*fco;
    }
  }
  SpectralDensity_sohmic::SpectralDensity_sohmic(){
    Parameters param;
    setup(param);
  }
  SpectralDensity_sohmic::SpectralDensity_sohmic(Parameters &param, const std::string &prefix){
    setup(param, prefix);
  }


///bump:
  std::string SpectralDensity_bump::add_name(const std::string &mpgname, const std::string &str)const{
    if(mpgname=="")return str;
    return std::string(mpgname+"_"+str);
  } 
  void SpectralDensity_bump::setup(Parameters &param, const std::string &prefix){

    double Emin = param.get_as_double(add_name(prefix,"E_min"),-hbar_in_meV_ps);
    double Emax = param.get_as_double(add_name(prefix,"E_max"), hbar_in_meV_ps);
    wmin = param.get_as_double(add_name(prefix,"omega_min"),Emin/hbar_in_meV_ps);
    wmax = param.get_as_double(add_name(prefix,"omega_max"),Emax/hbar_in_meV_ps);
 
    scale=param.get_as_double(add_name(prefix,"scale"),1.);
  }

  //NOTE: normalization to J(w=(wmax+wmin)/2) = 1/(2*pi), so that rate there is ~ 1
  double SpectralDensity_bump::f(double w)const{
    if(w<=wmin || w>=wmax)return 0.;
    double x=(w-wmin)/(wmax-wmin)*2.-1;
    return scale*exp(-x*x/(1.-x*x))/(2.*M_PI);
  }

  SpectralDensity_bump::SpectralDensity_bump(Parameters &param, const std::string &prefix){
    setup(param, prefix);
  }


}//namespace

