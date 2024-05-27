#include "ProcessTensorBuffer.hpp"
#include "Parameters.hpp"

using namespace ACE;

int main(int args, char **argv){
  Parameters param(args, argv, true);
  
  std::string read_PT=param.get_as_string_check("read_PT");
  int N_front=param.get_as_size_t("N_front");
  int N_back=param.get_as_size_t("N_back");
  
  if(N_front<2 && N_back<2){
    std::cerr<<"N_front<2 && N_back<2: nothing to extend"<<std::endl;
    exit(1);
  }

  ProcessTensorBuffer PTB(read_PT);
  PTB.expand_outer(N_front,N_back);

  return 0;
}
