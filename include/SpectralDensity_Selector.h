#ifndef ACE_SPECTRAL_DENSITY_SELECTOR_DEFINED_H
#define ACE_SPECTRAL_DENSITY_SELECTOR_DEFINED_H

#include "SpectralDensity.h"
#include "ReadTable.h"
#include "RealFunction_Interpolate.h"


class SpectralDensity_Selector{
public:
  RealFunctionPtr SD;

  operator RealFunctionPtr() { return SD; } 

  static std::string add_name(const std::string & prefix, const std::string & str){
    if(prefix=="")return str;
    return std::string(prefix+"_"+str);
  }
  static std::string get_prefix_without_J(const std::string & prefix=""){
    std::string prefix_without_J;
    if(prefix.size()>2 && 
       prefix[prefix.size()-2]=='_' &&  prefix[prefix.size()-1]=='J' ){
      prefix_without_J=prefix.substr(0, prefix.size()-2);
    }
    return prefix_without_J;
  }
  static bool is_specified(Parameters &param, const std::string & prefix=""){
    if(param.is_specified(add_name(prefix,"type")))return true;
    if(param.is_specified(add_name(prefix,"from_file")))return true;
    if(param.is_specified(get_prefix_without_J(prefix)+std::string("_rate")))return true;
    return false;
  }
  void setup(Parameters &param, const std::string & prefix=""){
    std::string rate_string=get_prefix_without_J(prefix)+std::string("_rate");

    RealFunctionPtr SD_; 
    std::string type=param.get_as_string(add_name(prefix,"type"));
    std::string fparam=add_name(prefix,"from_file");
    std::string file=param.get_as_string(fparam);
    if(file!=""){
      int col1=param.get_as_size_t(fparam, 0, 0, 1); //column containing omega
      int col2=param.get_as_size_t(fparam, 1, 0, 2); //column containing J
      SD_= new RealFunction_Interpolate(file, col1, col2);
    }else if(type=="" && !param.is_specified(rate_string)){
      std::cerr<<"Spectral density: Please specify either '"<<fparam<<"' or '"<<add_name(prefix,"type")<<"'!"<<std::endl;
      std::cerr<<"Parameters:"<<std::endl<<param;
      exit(1);

    }else if(type=="QDPhonon"){
      SD_ = new SpectralDensity_QD(param, prefix);

    }else if(type=="sohmic"||type=="ohmic"){
      SD_ = new SpectralDensity_sohmic(param, prefix);

    }else if(type=="const"){
      double val = param.get_as_double(add_name(prefix,"type"),0,0,1);
      SD_ = new RealFunction_Const(val);
     
    }else if(type=="flat_rate"||param.is_specified(rate_string) ){
      double rate;
      if(param.is_specified(rate_string)){
        rate = param.get_as_double(rate_string);
      }else{
        rate = param.get_as_double(add_name(prefix,"type"),0,0,1);
      }
      double j=rate/(2.*M_PI); 
      SD_ = new RealFunction_Const(j);

    }else if(type=="lorentzian"){
// sc and c can be set via .._J_scale and .._J_shift -> make more consistent:

      double gamma=param.get_as_double_check(add_name(prefix,"gamma"));
      SD_ = new RealFunction_Lorentzian(1., gamma, 0.);
    }else{
      std::cerr<<"Spectral density type '"<<type<<"' not known!"<<std::endl;
      exit(1);
    }
 
// Apply cut-offs 

    if(param.is_specified(add_name(prefix,"cutoff_max_logistic"))){
      if(param.get_nr_cols(add_name(prefix,"cutoff_max_logistic"),0)<2){
        std::cerr<<"Usage: cutoff_max_logistic  CENTER  WIDTH!"<<std::endl;
        exit(1);
      }
      double c=param.get_as_double(add_name(prefix,"cutoff_max_logistic"),0,0,0);
      double width=param.get_as_double(add_name(prefix,"cutoff_max_logistic"),1,0,1);
      RealFunctionPtr tmp=new RealFunction_Logistic(c, -width);
      RealFunctionPtr tmp2=new RealFunction_Product(SD_, tmp);
      SD_=tmp2;
    }
    if(param.is_specified(add_name(prefix,"cutoff_min_logistic"))){
      if(param.get_nr_cols(add_name(prefix,"cutoff_min_logistic"),0)<2){
        std::cerr<<"Usage: cutoff_min_logistic  CENTER  WIDTH!"<<std::endl;
        exit(1);
      }
      double c=param.get_as_double(add_name(prefix,"cutoff_min_logistic"),0,0,0);
      double width=param.get_as_double(add_name(prefix,"cutoff_min_logistic"),1,0,1);
      RealFunctionPtr tmp=new RealFunction_Logistic(c, width);
      RealFunctionPtr tmp2=new RealFunction_Product(SD_, tmp);
      SD_=tmp2;
    }
    if(param.is_specified(add_name(prefix,"cutoff_abs_logistic"))){
      if(param.get_nr_cols(add_name(prefix,"cutoff_abs_logistic"),0)<2){
        std::cerr<<"Usage: cutoff_abs_logistic  CENTER  WIDTH!"<<std::endl;
        exit(1);
      }
      double c=param.get_as_double(add_name(prefix,"cutoff_abs_logistic"),0,0,0);
      double width=param.get_as_double(add_name(prefix,"cutoff_abs_logistic"),1,0,1);
      RealFunctionPtr tmp1=new RealFunction_Logistic(c, -width);
      RealFunctionPtr tmp2=new RealFunction_Chain(tmp1, RealFunctionPtr_Abs);
      RealFunctionPtr tmp3=new RealFunction_Product(SD_, tmp2);
      SD_=tmp3;
    }


//------ Final shift and scale
    double scale=param.get_as_double(add_name(prefix,"scale"),1.);
    double shift=param.get_as_double(add_name(prefix,"shift"),0.);
    SD = new RealFunction_ScaleShift(SD_, scale, shift);

//------ Print if required
    std::string print_J=param.get_as_string(add_name(prefix,"print"));
    if(print_J!=""){
      std::vector<std::string> svec=param.get_row(add_name(prefix,"print"),0);
      double wa=0., we=10.; int Ndiscr=1000;
      if(svec.size()>1){wa=Reader::readDouble(svec[1],add_name(prefix,"print")+" ->wa<-  we  N " ); }
      if(svec.size()>2){we=Reader::readDouble(svec[2],add_name(prefix,"print")+" wa  ->we<-  N " ); }
      if(svec.size()>3){Ndiscr=Reader::readSizeT(svec[3],add_name(prefix,"print")+" wa we ->N<- " ); }
      SD->print(print_J, wa, we, Ndiscr);
    }
  }
  SpectralDensity_Selector(){
    Parameters param;
    setup(param);
  }
  SpectralDensity_Selector(Parameters &param, const std::string & prefix=""){
    setup(param, prefix);
  }
};

#endif
