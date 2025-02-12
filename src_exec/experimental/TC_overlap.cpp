#include "Trafo_Chain.hpp"
#include "Parameters.hpp"

using namespace ACE;

int main(int args, char ** argv){
  Parameters param(args, argv, false);
 
  std::vector<std::string> infile=param.get_all_strings("TC");
  if(infile.size()<2){
    std::cerr<<"'TC' needs two arguments!"<<std::endl;
    exit(1);
  }


 
  Trafo_Chain tc1(infile[0]);
  tc1.print_info();

  Trafo_Chain tc2(infile[1]);
  tc2.print_info();

  double compress_first=param.get_as_double("compress_first", 0);
  if(compress_first>0){
    std::cout<<"compressing TCs with epsilon="<<compress_first<<std::endl;
    tc1.compress(compress_first);
    tc1.print_info();
    tc2.compress(compress_first);
    tc2.print_info();
  }

  double overlap=tc1.overlap(tc2);
  std::cout<<"overlap: "<<overlap<<" of "<<tc1.lastdim()<<" or "<<tc2.lastdim()<<std::endl;

  return 0;
}
