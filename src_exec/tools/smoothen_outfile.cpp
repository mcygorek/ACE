#include "Constants.hpp"
#include "ReadTable.hpp"
#include "Parameters.hpp"
#include "ReaderBasics.hpp"
#include <fstream>

using namespace ACE;

int main(int args, char **argv){
  Parameters param(args, argv);

  std::string infile=param.get_as_string_check("infile");
  int col=param.get_as_size_t_check("col");
  double FWHM=param.get_as_double_check("FWHM");
  double n_sigma=param.get_as_double("n_sigma", 5.);
 
  bool do_append=param.get_as_bool("append", true);

  ReadTable tab(infile, 0, col-1);
  int N=tab.table.size();
  if(N<2){
    std::cerr<<"Need at least 2 data point!"<<std::endl;
    exit(1);
  }
  double ta=tab.table[0][0];
  double h=tab.table[1][0]-tab.table[0][0];
  if(tab.table[0][1] == 1./0.)tab.table[0][1]=0.;


  std::ifstream ifs(infile.c_str());

  std::string in_line;
  for(int i=0; i<N; i++){
    double res=0.;
    for(int j=0; j<N; j++){
      res+=h*gauss_from_FWHM(tab.table[j][0]-tab.table[i][0],FWHM)*tab.table[j][1];
    }

    if(!getRelevantLine(ifs, in_line)){
      std::cerr<<"Error: getRelevantLine(ifs, in_line) failed!"<<std::endl;
      exit(1);
    }
    
    if(do_append){ std::cout<<in_line; } else { std::cout<<ta+i*h; }
    std::cout<<" "<<res<<std::endl;

  }

  return 0;
}
