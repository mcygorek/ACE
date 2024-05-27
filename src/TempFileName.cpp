#include "TempFileName.hpp"
#include <string>
#include <fstream>
#include <iostream>
#include <random>
#include "DummyException.hpp"

namespace ACE{

void TempFileName::swap(TempFileName &other){
  std::string bck=fname;
  fname=other.fname;
  other.fname=bck;
}

std::string TempFileName::randomString(int length){
  static const char alphanum[] =
     "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_int_distribution<int> dist(0, sizeof(alphanum)-2); //note: sizeof also counts trailing '\0' and uniform_int_distribution uses range [start, end] with end inclusive

  std::string str(length,'\0');
  for(int i=0; i<length; i++){
    str[i]=alphanum[dist(mt)];
  }
  return str;
}

void TempFileName::initialize(const std::string &str){
//  fname = std::tmpnam(nullptr);
  noremove=false;
  int counter=0;
  bool exists;
  do{
    counter++;
    if(counter>=10){
      std::cerr<<"Can't create a new temporary file name!"<<std::endl;
      throw DummyException();
    }

    fname = "_tmp_"+str+"_"+randomString(8);
    std::ifstream f(fname.c_str());
    exists=f.good();

  }while(exists);
}
TempFileName & TempFileName::set_noremove(bool nr){
  noremove=nr;
  return *this;
}

void TempFileName::remove(){
  if(fname!="" && !noremove){
    std::remove(fname.c_str());
  }
}

}//namespace
