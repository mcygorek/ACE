#include "ReadExpression.hpp"
#include <iostream>

using namespace ACE;

int main(int args, char **argv){
  if(args<2){std::cerr<<"Need one argument!"<<std::endl;exit(1);}

  std::cout<<"input: '"<<argv[1]<<"'"<<std::endl;
  Eigen::MatrixXcd op=ReadExpression(argv[1],false);
//  Eigen::MatrixXcd op=ReadExpression(argv[1],true);
  std::cout<<"result:"<<std::endl<<op<<std::endl;

  return 0;
}

