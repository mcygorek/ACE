#include "Parameters.hpp"
#include "ProcessTensorBuffer.hpp"

using namespace ACE;

int main(int args, char** argv){
 try{ 
  Parameters param(args, argv);
  std::string read_PT=param.get_as_string_check("read_PT");

  ProcessTensorBuffer PTB(read_PT);
  for(int n=0; n<PTB.get_n_tot(); n++){
    ProcessTensorElement &e=PTB.get(n, ForwardPreload);
    e.env_ops.remove_all_but_first();
  }
 }catch(DummyException &e){  
  return 1;
 }
  return 0;
}
