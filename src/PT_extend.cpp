#include "InfluenceFunctional_OD.h"
#include "Parameters.h"


int main(int args, char **argv){
  Parameters param(args, argv, true);
  
  std::string read_PT=param.get_as_string("read_PT");
  std::string write_PT=param.get_as_string("write_PT");
  int N_front=param.get_as_size_t("N_front");
  int N_back=param.get_as_size_t("N_back");

  if(read_PT==""){
    std::cerr<<"Please specify 'read_PT'!"<<std::endl;
    exit(1);
  }
  if(write_PT==""){
    std::cerr<<"Please specify 'write_PT'!"<<std::endl;
    exit(1);
  }
  if(N_front<2 && N_back<2){
    std::cerr<<"N_front<2 && N_back<2: nothing to extend"<<std::endl;
    exit(1);
  }

  InfluenceFunctional_OD IF(read_PT);
  
  if(N_front>1)IF.dict.expand_space_front(N_front);
  if(N_back>1)IF.dict.expand_space_back(N_back);

  IF.write_binary(write_PT);
  
}
