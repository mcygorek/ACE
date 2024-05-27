#include "ReadTable.hpp"
#include "Parameters.hpp"
#include "ReaderBasics.hpp"
#include <fstream>

using namespace ACE;

int main(int args, char **argv){
  Parameters param(args, argv);

//  std::string outfile=param.get_as_string_check("outfile");
  std::string infile=param.get_as_string_check("infile");
  int col=param.get_as_size_t_check("col");
  bool halfstep=param.get_as_bool("halfstep", false);
 
  bool do_append=param.get_as_bool("append", true);

  bool use_Simpson=param.get_as_bool("use_Simpson",false);

  ReadTable tab(infile, 0, col-1);
  int N=tab.table.size();


 if(use_Simpson){
  if(N<3){
    std::cerr<<"Simpson integration requires at least 3 data points!"<<std::endl;
    exit(1);
  }

  double ta=tab.table[0][0];
  double h=tab.table[1][0]-tab.table[0][0];
//  std::cerr<<"h: "<<h<<std::endl;
  if(tab.table[0][1] == 1./0.)tab.table[0][1]=0.;
  double last=param.get_as_double("init",0.);

  int N0=(N-1)/2;
  std::ifstream ifs;
  std::string in_line;
  if(do_append)ifs.open(infile.c_str());
 
  if(do_append && !getRelevantLine(ifs, in_line)){
    std::cerr<<"Error: getRelevantLine(ifs, in_line) failed!"<<std::endl;
    exit(1);
  }
  if(do_append){ std::cout<<in_line; } else { std::cout<<ta; }
  std::cout<<" "<<last<<std::endl;

  for(int i=0; i<N0; i++){
    if(do_append && !getRelevantLine(ifs, in_line)){
      std::cerr<<"Error: getRelevantLine(ifs, in_line) failed!"<<std::endl;
      exit(1);
    }

    if(halfstep){
      double res = last + h/12.*(  \
      -tab.table[2*i][1] + 8.*tab.table[2*i+1][1]  + 5.*tab.table[2*i+2][1]);

      if(do_append){ std::cout<<in_line; } else { std::cout<<ta+(2.*i+1)*h; }
      std::cout<<" "<<res<<std::endl;
    }


    double res = last +  h/3.*(    \
          tab.table[2*i][1] + 4.*tab.table[2*i+1][1]  + tab.table[2*i+2][1]);

    if(!getRelevantLine(ifs, in_line)){
      std::cerr<<"Error: getRelevantLine(ifs, in_line) failed!"<<std::endl;
      exit(1);
    }
    
    if(do_append){ std::cout<<in_line; } else { std::cout<<ta+(2.*i+1)*h; }
    std::cout<<" "<<res<<std::endl;

    last=res;
  }

  if(N>2*N0+1){
    double res=last+h*tab.table[2*N0+1][1];
    std::cout<<ta+(2*N0+2)*h<<" "<<res<<std::endl;
  }


 }else{  //use simple midpoint rule

  if(N<2){
    std::cerr<<"integration requires at least 2 data points!"<<std::endl;
    exit(1);
  }

  std::ifstream ifs;
  if(do_append)ifs.open(infile.c_str());

  double last=param.get_as_double("init",0.);

  for(int n=0; n<N; n++){
    std::string in_line;
    if(do_append && !getRelevantLine(ifs, in_line)){
      std::cerr<<"Error: getRelevantLine(ifs, in_line) failed!"<<std::endl;
      exit(1);
    }

    if(n>0){
      double h=tab.table[n][0]-tab.table[n-1][0];
      last+=0.5*h*(tab.table[n][1]+tab.table[n-1][1]);
    }
    
    if(do_append){ std::cout<<in_line; } else { std::cout<<tab.table[n][0]; }
    std::cout<<" "<<last<<std::endl;
  }
 }

  return 0;
}
