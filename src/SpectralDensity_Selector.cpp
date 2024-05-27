#include "SpectralDensity_Selector.hpp"
#include "ReadTable.hpp"
#include "RealFunction_Interpolate.hpp"
#include "Parameters.hpp"
#include "Reader.hpp"
#include "DummyException.hpp"

namespace ACE{

  std::string SpectralDensity_Selector::add_name(const std::string & prefix, const std::string & str){
    if(prefix=="")return str;
    return std::string(prefix+"_"+str);
  }

  std::string SpectralDensity_Selector::get_prefix_without_J(const std::string & prefix){
    std::string prefix_without_J;
    if(prefix.size()>2 && 
       prefix[prefix.size()-2]=='_' &&  prefix[prefix.size()-1]=='J' ){
      prefix_without_J=prefix.substr(0, prefix.size()-2);
    }
    return prefix_without_J;
  }

  bool SpectralDensity_Selector::is_specified(Parameters &param, const std::string & prefix){
    if(param.is_specified(add_name(prefix,"type")))return true;
    if(param.is_specified(add_name(prefix,"from_file")))return true;
    if(param.is_specified(get_prefix_without_J(prefix)+std::string("_rate")))return true;
    return false;
  }

  void SpectralDensity_Selector::setup(Parameters &param, const std::string & prefix){
    std::string rate_string=get_prefix_without_J(prefix)+std::string("_rate");

    RealFunctionPtr SD_; 
    std::string type=param.get_as_string(add_name(prefix,"type"));
    std::string fparam=add_name(prefix,"from_file");
    std::string file=param.get_as_string(fparam);
    if(file!=""){
      std::cout<<"Using spectral density 'from_file'"<<std::endl<<std::endl;
      int col1=param.get_as_size_t(fparam, 0, 0, 1); //column containing omega
      int col2=param.get_as_size_t(fparam, 1, 0, 2); //column containing J
      SD_= std::make_shared<RealFunction_Interpolate>(file, col1, col2);
    }else if(type=="" && !param.is_specified(rate_string)){
      std::cerr<<"Spectral density: Please specify either '"<<fparam<<"' or '"<<add_name(prefix,"type")<<"'!"<<std::endl<<std::endl;
      throw DummyException();

    }else if(type=="QDPhonon"){
      std::cout<<"Using spectral density 'QDPhonon'"<<std::endl;
      SD_ = std::make_shared<SpectralDensity_QD>(param, prefix);

    }else if(type=="sohmic"||type=="ohmic"){
      std::cout<<"Using spectral density '(s-)ohmic'"<<std::endl;
      SD_ = std::make_shared<SpectralDensity_sohmic>(param, prefix);

    }else if(type=="bump"||type=="Bump"){
      std::cout<<"Using spectral density 'bump'"<<std::endl;
      SD_ = std::make_shared<SpectralDensity_bump>(param, prefix);

    }else if(type=="const"){
      std::cout<<"Using spectral density 'const'"<<std::endl;
      double val = param.get_as_double(add_name(prefix,"type"),0,0,1);
      SD_ = std::make_shared<RealFunction_Const>(val);
     
    }else if(type=="flat_rate"||param.is_specified(rate_string) ){
      std::cout<<"Using spectral density 'rate'"<<std::endl;
      double rate;
      if(param.is_specified(rate_string)){
        rate = param.get_as_double(rate_string);
      }else{
        rate = param.get_as_double(add_name(prefix,"type"),0,0,1);
      }
      double j=rate/(2.*M_PI); 
      SD_ = std::make_shared<RealFunction_Const>(j);

    }else if(type=="lorentzian"){
      std::cout<<"Using spectral density 'lorentzian'"<<std::endl;
// sc and c can be set via .._J_scale and .._J_shift -> make more consistent:

      double gamma=param.get_as_double_check(add_name(prefix,"gamma"));
      SD_ = std::make_shared<RealFunction_Lorentzian>(1., gamma, 0.);
    }else{
      std::cerr<<"Spectral density type '"<<type<<"' not known!"<<std::endl;
      throw DummyException();
    }
 
// Apply cut-offs 

    if(param.is_specified(add_name(prefix,"cutoff_max_logistic"))){
      if(param.get_nr_cols(add_name(prefix,"cutoff_max_logistic"),0)<2){
        std::cerr<<"Usage: cutoff_max_logistic  CENTER  WIDTH!"<<std::endl;
        throw DummyException();
      }
      double c=param.get_as_double(add_name(prefix,"cutoff_max_logistic"),0,0,0);
      double width=param.get_as_double(add_name(prefix,"cutoff_max_logistic"),1,0,1);
      if(fabs(width)<1e-36){
        std::cerr<<"Please set cutoff_max_logistic width !=0"<<std::endl;
        throw DummyException();
      }
      RealFunctionPtr tmp=std::make_shared<RealFunction_Logistic>(c, -width);
      RealFunctionPtr tmp2=std::make_shared<RealFunction_Product>(SD_, tmp);
      SD_=tmp2;
    }
    if(param.is_specified(add_name(prefix,"cutoff_min_logistic"))){
      if(param.get_nr_cols(add_name(prefix,"cutoff_min_logistic"),0)<2){
        std::cerr<<"Usage: cutoff_min_logistic  CENTER  WIDTH!"<<std::endl;
        throw DummyException();
      }
      double c=param.get_as_double(add_name(prefix,"cutoff_min_logistic"),0,0,0);
      double width=param.get_as_double(add_name(prefix,"cutoff_min_logistic"),1,0,1);
      if(fabs(width)<1e-36){
        std::cerr<<"Please set cutoff_min_logistic width !=0"<<std::endl;
        throw DummyException();
      }
      RealFunctionPtr tmp=std::make_shared<RealFunction_Logistic>(c, width);
      RealFunctionPtr tmp2=std::make_shared<RealFunction_Product>(SD_, tmp);
      SD_=tmp2;
    }
    if(param.is_specified(add_name(prefix,"cutoff_abs_logistic"))){
      if(param.get_nr_cols(add_name(prefix,"cutoff_abs_logistic"),0)<2){
        std::cerr<<"Usage: cutoff_abs_logistic  CENTER  WIDTH!"<<std::endl;
        throw DummyException();
      }
      double c=param.get_as_double(add_name(prefix,"cutoff_abs_logistic"),0,0,0);
      double width=param.get_as_double(add_name(prefix,"cutoff_abs_logistic"),1,0,1);
      if(fabs(width)<1e-36){
        std::cerr<<"Please set cutoff_abs_logistic width !=0"<<std::endl;
        throw DummyException();
      }
      RealFunctionPtr tmp1=std::make_shared<RealFunction_Logistic>(c, -width);
      RealFunctionPtr tmp2=std::make_shared<RealFunction_Chain>(tmp1, RealFunctionPtr_Abs);
      RealFunctionPtr tmp3=std::make_shared<RealFunction_Product>(SD_, tmp2);
      SD_=tmp3;
    }

//------ Final shift and scale
    double scale=param.get_as_double(add_name(prefix,"scale"),1.);
    double shift=param.get_as_double(add_name(prefix,"shift"),0.);
    SD = std::make_shared<RealFunction_ScaleShift>(SD_, scale, shift);

//------ Print if required
    //use either of, e.g., 'Boson_J_print' or 'Boson_print_J'
    std::string J_print=param.get_as_string(add_name(prefix,"print"));
    std::string print_J=param.get_as_string(add_name(get_prefix_without_J(prefix),"print_J"));
    std::string J_print_par=add_name(get_prefix_without_J(prefix),"print_J");
    if(J_print!=""){
      J_print_par=add_name(prefix,"print");
      print_J=J_print;
    }

    if(print_J!=""){
      std::vector<std::string> svec=param.get_row(J_print_par,0);
      double wa=0., we=10.; int Ndiscr=1000;
      if(svec.size()>1){wa=readDouble(svec[1],J_print_par+" ->wa<-  we  N " ); }
      if(svec.size()>2){we=readDouble(svec[2],J_print_par+" wa  ->we<-  N " ); }
      if(svec.size()>3){Ndiscr=readSizeT(svec[3],J_print_par+" wa we ->N<- " ); }
      SD->print(print_J, wa, we, Ndiscr);

      std::cout<<"Printing spectral density to file '"<<print_J<<"'"<<std::endl;
    }
  }
  SpectralDensity_Selector::SpectralDensity_Selector(){
    Parameters param;
    setup(param);
  }

}//namespace
