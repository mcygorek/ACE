#include "Parameters.hpp"
#include "ReadPT_struct.hpp"

namespace ACE{

void ReadPT_struct::setup(Parameters &param, const std::string &par_name, const std::string & expand_name){
  if(par_name=="" || expand_name==""){ 
    std::cerr<<"ReadPT_struct::setup: par_name==\"\" || expand_name==\"\" !"<<std::endl;
    exit(1);
  }
  expand_front=1; 
  expand_back=1;
  if(param.is_specified(par_name)){
    if(param.is_specified(expand_name)){
      std::cerr<<"ReadPT_struct::setup: Please specify only one of '"<<par_name<<"' and '"<<expand_name<<"'!"<<std::endl;
      exit(1);
    }
    fname=param.get_as_string_check(par_name);
  }else if(param.is_specified(expand_name)){
    param.complain_if_row_shorter(expand_name, 3, 0, expand_name+" FILENAME expand_front expand_back");
    fname=param.get_as_string_check(expand_name);
    expand_front=param.get_as_double_check(expand_name, 0, 1); 
    expand_back=param.get_as_double_check(expand_name, 0, 2); 
  }else{
    fname="";
  }
}

void ReadPT_struct::print_info(std::ostream &ofs)const{
  ofs<<"fname="<<fname<<" expand_front="<<expand_front<<" expand_back="<<expand_back<<std::endl;
}
}//namespace
