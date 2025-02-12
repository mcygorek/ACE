#include "Pulse_Selector.hpp"
#include "Operators.hpp"
#include "Pulse.hpp"
#include "ReadExpression.hpp"
#include "Reader.hpp"
#include "Constants.hpp"
#include "ComplexFunction_Interpolate.hpp"
#include "DummyException.hpp"
#include <iostream>

namespace ACE{

std::pair<ComplexFunctionPtr, Eigen::MatrixXcd>  Pulse_Selector(
                                      const std::vector<std::string> &toks){

  if(toks.size()<1){
    std::cerr<<"Pulse_Selector: toks.size()<1!"<<std::endl;
    throw DummyException();
  }

  std::string type=toks[0];

  if(type=="Gauss"||type=="Gaussian"){ //center, FWHM, area, detuning
    std::string errmsg="Usage: add_Pulse Gauss CENTER FWHM AREA DETUNING OPERATOR!";
    if(toks.size()<4){
      std::cerr<<"Not enough parameters: "<<errmsg<<std::endl;
      throw DummyException();
    }
    //defaults:
    double detuning=(toks.size()<5) ? 0 
           : readDouble(toks[4],errmsg+"\nPulse: Gauss: detuning");

    Eigen::MatrixXcd Op=(toks.size()<6) 
           ?  hbar_in_meV_ps/2.*Operators(2).ketbra(1,0)
           :  (Eigen::MatrixXcd) ReadExpression(toks[5]);

       
    ComplexFunctionPtr pulse=std::make_shared<Pulse_Gauss>(
      readDouble(toks[1],errmsg+"\nPulse: Gauss: center(time)"),
      readDouble(toks[2],errmsg+"\nPulse: Gauss: FWHM(time)"),
      readDouble(toks[3],errmsg+"\nPulse: Gauss: area"),
      detuning);
   
    return std::make_pair(pulse, Op);

  }else if(type=="rect"){
 
    std::string errmsg="Usage: add_Pulse rect T_ON T_OFF T_RISE[=0] SCALE[=1] DETUNING[=0] OPERATOR[={hbar/2*|1><0|_2}]!";
    if(toks.size()<3){
      std::cerr<<"Not enough parameters: "<<errmsg<<std::endl;
      throw DummyException();
    }

    double t_rise=(toks.size()<4) ? 0 
           : readDouble(toks[3],errmsg+"\nPulse: rect: t_rise(ps)");
    
    double scale=(toks.size()<5) ? 1. 
           : readDouble(toks[4],errmsg+"\nPulse: rect: scale");
          
    double detuning=(toks.size()<6) ? 0
           : readDouble(toks[5],errmsg+"\nPulse: rect: detuning(meV)");

    Eigen::MatrixXcd Op=(toks.size()<7)
           ?  hbar_in_meV_ps/2.*Operators(2).ketbra(1,0)
           :  (Eigen::MatrixXcd) ReadExpression(toks[6]);
    

    ComplexFunctionPtr pulse=std::make_shared<Pulse_SmoothRect>(
           readDouble(toks[1],errmsg+"\nPulse: rect: t_on(ps)"), 
           readDouble(toks[2],errmsg+"\nPulse: rect: t_off(ps)"), 
           t_rise, scale, detuning);

   return std::make_pair(pulse, Op);

  }else if(type=="file"||type=="File"){
    std::string errmsg="Usage: add_Pulse File FILENAME [SYSOP={hbar/2*|1><0|_2}] !";
    if(toks.size()<2){
      std::cerr<<"Not enough parameters: "<<errmsg<<std::endl;
      throw DummyException();
    }
    std::string fname=toks[1];
    
    std::string OpString="{hbar/2*|1><0|_2}";
    if(toks.size()>2){
      OpString=toks[2];
    }
          
    int tcol=0, col=1, im_col=col+1;

    Eigen::MatrixXcd expr;
    try{
      ComplexFunctionPtr pulse =
       std::make_shared<ComplexFunction_Interpolate>(fname, tcol, col, im_col); 
      return std::make_pair(pulse, ReadExpression(OpString));
    }catch(DummyException &e){
      std::cerr<<"called by: add_Pulse "<<fname<<" "<<OpString<<std::endl;
      throw e;
    }
  }else{
    std::cerr<<"Cannot recognize pulse type '"<<type<<"'!"<<std::endl;
    throw DummyException();
  }
}

}//namespace
