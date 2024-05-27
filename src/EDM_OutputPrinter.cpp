#include "EDM_OutputPrinter.hpp"
#include <iomanip>

namespace ACE{

void EDM_OutputPrinter::print_time(double t){
  if(ofs){
    (*ofs)<<t;
  }
}
void EDM_OutputPrinter::print_value(std::complex<double> val){
  if(ofs){
    (*ofs)<<" "<<val.real()<<" "<<val.imag();
  }
}
void EDM_OutputPrinter::print_endl(){
  if(ofs){
    (*ofs)<<std::endl;
  }
}

void EDM_OutputPrinter::set_file(const std::string &str){
  if(ofs){
    if(ofs->is_open())ofs->close();
    ofs.reset(nullptr);
  }
  if(str!="" && str!="/dev/null"){
    ofs.reset(new std::ofstream(str.c_str()));
  }
}
void EDM_OutputPrinter::setup(Parameters &param, const std::string &outfile_default, const std::string &outfile_key){

  std::string outfile=param.get_as_string(outfile_key, outfile_default);
  set_file(outfile);

  int set_precision=param.get_as_int("set_precision",-1);
  if(ofs && set_precision>0)*ofs<<std::setprecision(set_precision);
}

}//namespace
