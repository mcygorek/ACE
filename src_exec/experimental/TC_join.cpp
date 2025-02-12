#include "ACE.hpp"

using namespace ACE;

int main(int args, char ** argv){
  Parameters param(args, argv, true);
 
  std::vector<std::string> infile=param.get_all_strings("TC");
  if(infile.size()<1){
    std::cerr<<"'TC' needs at least one arguments!"<<std::endl;
    exit(1);
  }

  double print_ortho_threshold=param.get_as_double("print_ortho_threshold",1e-8);
  std::string outfile=param.get_as_string_check("outfile");
  
  Trafo_Chain tc(infile[0]);
  std::cout<<"TC[0]: ";
  tc.print_info();
  if(tc.size()<1){
    std::cerr<<"TC[0].size()<1!"<<std::endl;
    exit(1);
  }
   
  double compress=param.get_as_double("compress",0);
  double compress_first=param.get_as_double("compress_first",0);
  if(compress_first>0){
    std::cout<<"compressing first TC with epsilon="<<compress_first<<std::endl;
    tc.compress(compress_first);
  }

  
  double threshold=param.get_as_double("threshold");
  int max_add=param.get_as_int("max_add",-1);
  for(size_t i=1; i<infile.size(); i++){
    Trafo_Chain tc2(infile[i]);
    if(compress_first>0)tc2.compress(compress_first);
    std::cout<<"TC["<<i<<"]: ";
    tc2.print_info();
    {Eigen::MatrixXd OM=tc2.overlap_matrix(tc2);
    std::cout<<"Overlaps: max_diff_from_ortho="<<max_diff_from_ortho(OM)<<std::endl;}
    std::cout<<"old overlap: "<<tc.overlap(tc2)<<" of "<<tc.lastdim()<<" or "<<tc2.lastdim()<<std::endl;

    tc.add_ortho(tc2, threshold, compress, max_add);

    {Eigen::MatrixXd OM=tc.overlap_matrix(tc);
    std::cout<<"Overlaps: max_diff_from_ortho="<<max_diff_from_ortho(OM)<<std::endl;
    print_diff_from_ortho(OM, print_ortho_threshold);}



    std::cout<<"new overlap: "<<tc.overlap(tc2)<<" of "<<tc.lastdim()<<" or "<<tc2.lastdim()<<std::endl;
    std::cout<<"TC[0]: ";
    tc.print_info();
  }

  tc.write(outfile);

  return 0;
}
