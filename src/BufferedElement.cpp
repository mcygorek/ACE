#include "BufferedElement.hpp"
#include "DummyException.hpp"
#include <sstream>
#include <iostream>

namespace ACE{


BufferedElement::~BufferedElement(){};

void BufferedInt::read_binary(std::istream &is){
  std::string line;
  std::getline(is, line);
  std::istringstream iss(line);
  if(!(iss>>value)){
    std::cerr<<"Can't read BufferedInt!"<<std::endl;
    throw DummyException();
  }
}
void BufferedInt::write_binary(std::ostream &os)const{
  os<<value<<" "<<std::endl;
} 
BufferedInt::~BufferedInt(){};

}//namespace
