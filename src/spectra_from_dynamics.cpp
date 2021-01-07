#include "ReadTable.h"
#include "Parameters.h"
#define SLOWFT_PRINT
#include "slowFT.h"
#include "FT_Output.h"

int main(int args, char** argv){
 
  Parameters param(args, argv, true); 
 
  FT_Output fto(param);

  std::string dynamics_outfile=param.get_as_string("outfile");
  if(dynamics_outfile==""){
    std::cerr<<"Please specify input (dynamics) as 'outfile'!"<<std::endl;
    exit(1);
  }
  Simulation_Results res(dynamics_outfile);
 
//  res.print("TEST.res");

  std::cout<<"res.size(): "<<res.size()<<std::endl;
  if(res.size()<2){std::cerr<<"res.size()<3!"<<std::endl; exit(1);}
  std::cout<<"res[0].second.size(): "<<res[0].second.size()<<std::endl;
  
  fto.print(res);

  return 0;
}
