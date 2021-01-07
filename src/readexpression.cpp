#include "ReadExpression.h"



int main(int args, char **argv){
  if(args<2){std::cerr<<"Need one argument!"<<std::endl;exit(1);}

  std::cout<<"input: '"<<argv[1]<<"'"<<std::endl;
  std::complex<double> c=ReadExpression(argv[1]);
  std::cout<<"result: "<<c.real()<<" "<<c.imag()<<std::endl;


  return 0;
}
