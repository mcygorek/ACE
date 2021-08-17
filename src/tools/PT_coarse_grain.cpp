#include "Modify_PT.h"
#include "Parameters.h"

int main(int args, char **argv){
  Parameters param(args, argv, true, false);
  
  std::string read_PT=param.get_as_string_check("read_PT");
  std::string write_PT=param.get_as_string_check("write_PT");
  int n_coarse=param.get_as_size_t_check("steps");
  double dict_zero=param.get_as_double("dict_zero", 1e-20);

  InfluenceFunctional_OD IF(read_PT);

  Modify_PT::coarse_grain(IF, n_coarse, dict_zero);

  IF.write_binary(write_PT);
  
  std::cout<<"PT written to file '"<<write_PT<<"'!"<<std::endl;
  exit(1);
}
