#include "Modify_PT.h"
#include "Parameters.h"

int main(int args, char **argv){
  Parameters param(args, argv, true, true);
  
  double ta=param.get_as_double("ta", 0);
  double dt=param.get_as_double("dt", 1e-2);

  std::string read_PT=param.get_as_string_check("read_PT");
  std::string write_PT=param.get_as_string_check("write_PT");

  int n_coarse=param.get_as_size_t("n_coarse");
  double dict_zero=param.get_as_double("dict_zero", 1e-20);

  InfluenceFunctional_OD IF(read_PT);
  FreePropagator prop(param);

  Modify_PT::apply_system_propagator(IF, prop, ta, dt, dict_zero);

  if(n_coarse>1){
    Modify_PT::coarse_grain(IF, n_coarse, dict_zero);
  }

  IF.write_binary(write_PT);
  
  std::cout<<"PT written to file '"<<write_PT<<"'!"<<std::endl;
  exit(1);
}
