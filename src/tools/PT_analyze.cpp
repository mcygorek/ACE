#include "Parameters.h"
#include "AnalyzePT.h"

int main(int args, char ** argv){
  Parameters param(args, argv, true, false);

  std::string outfile=param.get_as_string("outfile","PT_analyse.out");

  std::string read_PT=param.get_as_string("read_PT");
  if(read_PT==""){ 
    std::cerr<<"Please specify parameter 'read_PT'!"<<std::endl;
    exit(1);
  }

  InfluenceFunctional_OD IF(read_PT);

  AnalyzePT::print_summary(outfile, IF);
  
  std::string full_pgm=param.get_as_string("full_pgm");
  if(full_pgm!="")AnalyzePT::print_pgm(full_pgm, IF);
  
  std::vector<std::vector<std::string> > svv=param.get("single_pgm");
  for(int r=0; r<(int)svv.size(); r++){
    if(svv[r].size()<3){
      std::cerr<<"Usage: single_pgm FILENAME  n  beta!"<<std::endl;
      exit(1);
    }
    std::string fname=svv[r][0];
    int n=Reader::readSizeT(svv[r][1],"single_pgm: n");
    int a=Reader::readSizeT(svv[r][2],"single_pgm: beta");
    AnalyzePT::print_single_pgm(fname, IF, n, a);
  }

  
  AnalyzePT::print_summary(std::cout, IF);

  return 0;
}
