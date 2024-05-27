#include "ParametersScan.hpp"
#include "Parameters.hpp"

namespace ACE{

Parameters ParametersScan(Parameters &param, int scan, int verbosity){
  Parameters param2=param;
  std::string scan_name="scan"+std::to_string(scan);
  param2.add_from_prefix(scan_name, param);


  if(verbosity){Parameters param3; 
  param3.add_from_prefix(scan_name, param); 
  std::cout<<"added parameters for this scan: "<<param3.map.size()<<":"<<std::endl; 
  param3.print();std::cout<<std::endl;}


  std::string outfile=param.get_as_string("outfile", "ACE.out");
  std::string this_outfile=outfile;
  if(param.is_specified(scan_name+"_outfile")){
    this_outfile=param.get_as_string(scan_name+"_outfile");
  }else{
    this_outfile=outfile+"_"+scan_name;
  }
  param2.override_param("outfile", this_outfile);
  return param2;
}


}//namespace
