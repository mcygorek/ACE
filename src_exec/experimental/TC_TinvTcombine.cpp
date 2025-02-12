#include "ACE.hpp"

using namespace ACE;

int main(int args, char ** argv){
  Parameters param(args, argv, true);
 
  std::string outfile=param.get_as_string_check("outfile");
  std::vector<std::string> infile=param.get_all_strings("TC");
  if(infile.size()<1){
    std::cerr<<"'TC' needs at least one arguments!"<<std::endl;
    exit(1);
  }

  Trafo_Chain tc(infile[0]);
  std::cout<<"TC[0]: ";
  tc.print_info();
  if(tc.size()<1){
    std::cerr<<"TC[0].size()<1!"<<std::endl;
    exit(1);
  }
 
  double epsilon=param.get_as_double_check("epsilon");
 

  for(size_t i=1; i<infile.size(); i++){
    Trafo_Chain tc2(infile[i]);
    std::cout<<"TC["<<i<<"]: ";
    tc2.print_info();
    tc.combine(tc2);
  }

  tc.combine_Tinv_T(epsilon);

  tc.write(outfile);

  return 0;
}
