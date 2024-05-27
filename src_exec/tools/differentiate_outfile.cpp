#include "ReadTable.hpp"
#include "Parameters.hpp"
#include "ReaderBasics.hpp"
#include <fstream>
#include "DummyException.hpp"

using namespace ACE;

int main(int args, char **argv){
 try{
  Parameters param(args, argv);

  std::string infile=param.get_as_string_check("infile");
  std::string outfile=param.get_as_string("outfile");
  int col=param.get_as_size_t_check("col");
 
  bool do_append=param.get_as_bool("append", true);

  ReadTable tab(infile, 0, col-1);
  int N=tab.table.size();
  if(N<3){
    std::cerr<<"Central differences needs at least 3 data point!"<<std::endl;
    throw DummyException();
  }
  double ta=tab.table[0][0];
  double h=tab.table[1][0]-tab.table[0][0];
//  std::cerr<<"h: "<<h<<std::endl;
  if(tab.table[0][1] == 1./0.)tab.table[0][1]=0.;


  std::ifstream ifs(infile.c_str());

  std::string in_line;
  if(!getRelevantLine(ifs, in_line)){
    std::cerr<<"Error: getRelevantLine(ifs, in_line) failed!"<<std::endl;
    throw DummyException();
  }

  std::ostream *os; 
  std::ofstream filestream;
  if(outfile!="" && outfile!="/dev/null"){
    filestream.open(outfile.c_str());
    os=&filestream;
  }else{
    os=&std::cout;
  }
  

  for(int i=1; i<N-1; i++){

    double res = (tab.table[i+1][1] - tab.table[i-1][1])/(2.*h);

    if(!getRelevantLine(ifs, in_line)){
      std::cerr<<"Error: getRelevantLine(ifs, in_line) failed!"<<std::endl;
      throw DummyException();
    }
    
    if(do_append){ (*os)<<in_line; } else { (*os)<<ta+i*h; }
    (*os)<<" "<<res<<std::endl;

  }
  return 0;
 }catch(DummyException &e){
  return 1;
 }
}
