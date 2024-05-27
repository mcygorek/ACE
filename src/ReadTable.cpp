#include "ReadTable.hpp"
#include "Reader.hpp"
#include "DummyException.hpp"
#include <iostream>

namespace ACE{

  std::vector<double> & ReadTable::operator[] (int i){
    if(i>=(int)table.size()){
      std::cerr<<"ReadTable::operator[] out of bounds!"<<std::endl;
      throw DummyException();
    }
    return table[i];
  }
  
  int ReadTable::colums_in_first_row(const std::string &fname){
    std::ifstream ifs(fname.c_str());
    std::string firstline;
    if(!getRelevantLine(ifs,firstline)){
      std::cerr<<"Cannot read first line of file '"<<fname<<"'!"<<std::endl; 
      throw DummyException();
    }
    return number_of_doubles(firstline);
  }

  void ReadTable::read(const std::string &str, std::vector<int> indices){
    table.clear();
    std::ifstream is(str.c_str());
    if(!is.is_open()){ 
      std::cerr<<"ReadTable::read: cannot open file '"<<str<<"'!"<<std::endl;
      throw DummyException();
    }
    std::string line;
    while(std::getline(is, line)){
      std::vector<std::string> toks;
      tokenize(line, toks);
      if(toks.size()<1)continue;
      for(size_t i=0; i<indices.size(); i++){
        if((int)toks.size()<indices[i]+1){
          std::cerr<<"ReadTable::read: toks.size()<col["<<indices[i]<<"]!"<<std::endl;
          throw DummyException();
        }
      }
      std::vector<double> dv(indices.size());
      for(size_t i=0; i<indices.size(); i++){
        std::stringstream ss;
        ss<<"File: "<<str<<" column: "<<indices[i];
        dv[i]=readDouble(toks[indices[i]], ss.str());
      }
      table.push_back(dv);
    }
  }
  
  void ReadTable::read(const std::string &str, int col1, int col2, int col3){
    std::vector<int> vi(3);
    vi[0]=col1; 
    vi[1]=col2;
    vi[2]=col3;
    read(str,vi);
  }
  void ReadTable::read(const std::string &str, int col1, int col2){
    std::vector<int> vi(2);
    vi[0]=col1; 
    vi[1]=col2;
    read(str,vi);
  }
  void ReadTable::read(const std::string &str, int col){
    std::vector<int> vi(1,col);
    read(str,vi);
  }
  void ReadTable::read(const std::string &str){
    int ncol=colums_in_first_row(str);
    std::vector<int> indices;
    for(int i=0; i<ncol; i++)indices.push_back(i);
    read(str, indices);
  }


}//namespace
