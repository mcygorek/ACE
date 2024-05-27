#include "Parameters.hpp"
#include "Reader.hpp"
#include "ProcessTensorBuffer.hpp"
#include "DummyException.hpp"

using namespace ACE;

int main(int args, char ** argv){
//  Parameters param(args, argv, true, false);
  Parameters param(args, argv);

  std::string read_PT=param.get_as_string_check("read_PT");
  std::string write_PT=param.get_as_string_check("write_PT");

//  bool RWA=param.get_as_bool("RWA",false);

  ProcessTensorBuffer PTB(read_PT, true);

  int coarse_grain=param.get_as_size_t_check("coarse_grain");
  int buffer_blocksize=param.get_as_int("buffer_blocksize",PTB.blocksize);

  ProcessTensorBuffer PTB2;
  PTB2.set_new_file(write_PT,buffer_blocksize);
  PTB2.set_from_coarse_grain(PTB, coarse_grain); 

}
