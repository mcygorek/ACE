#include "Parameters.hpp"
#include "ProcessTensorBuffer.hpp"
#include "IF_OD_Dictionary.hpp"

using namespace ACE;

int main(int args, char** argv){
 try{ 
  Parameters param(args, argv);
//  std::string outfile=param.get_as_string_check("outfile");

  param.complain_if_row_shorter("read_PT_expand", 3, 0, 
                      "Usage: -read_PT_expand FILENAME dim_front dim_back!");

  ReadPT_struct readPT(param, "read_PT", "read_PT_expand");
  readPT.print_info();
  
  ProcessTensorBuffer PTB(readPT.fname);
  for(int n=0; n<PTB.get_n_tot(); n++){
    ProcessTensorElement &e=PTB.get(n, ForwardPreload);
    IF_OD_Dictionary & dict=e.accessor.dict;
    dict.expand(readPT,false);
  }
 }catch(DummyException &e){
  return 1;
 }
  return 0;
}
