#ifndef READ_TABLE_DEFINED_H
#define READ_TABLE_DEFINED_H

#include "Reader.h"

class ReadTable{
public:
  std::vector<std::vector<double> > table;

  size_t size()const{
    return table.size();
  }
  std::vector<double> & operator[] (int i){
    if(i>=(int)table.size()){
      std::cerr<<"ReadTable::operator[] out of bounds!"<<std::endl;
      exit(1);
    }
    return table[i];
  }
  operator std::vector<std::vector<double> > () const{
    return table;
  }
  static int colums_in_first_row(const std::string &fname){
    std::ifstream ifs(fname.c_str());
    std::string firstline;
    if(!Reader::getRelevantLine(ifs,firstline)){
      std::cerr<<"Cannot read first line of file '"<<fname<<"'!"<<std::endl; 
      exit(1);
    }
    return Reader::number_of_doubles(firstline);
  }
  void read(const std::string &str, std::vector<int> indices){
    table.clear();
    std::ifstream is(str.c_str());
    std::string line;
    while(std::getline(is, line)){
      std::vector<std::string> toks;
      Reader::tokenize(line, toks);
      if(toks.size()<1)continue;
      for(size_t i=0; i<indices.size(); i++){
        if((int)toks.size()<indices[i]+1){
          std::cerr<<"ReadTable::read: toks.size()<col["<<indices[i]<<"]!"<<std::endl;
          exit(1);
        }
      }
      std::vector<double> dv(indices.size());
      for(size_t i=0; i<indices.size(); i++){
        std::stringstream ss;
        ss<<"File: "<<str<<" column: "<<indices[i];
        dv[i]=Reader::readDouble(toks[indices[i]], ss.str());
      }
      table.push_back(dv);
    }
  }  
  void read(const std::string &str, int col1, int col2, int col3){
    std::vector<int> vi(3);
    vi[0]=col1; 
    vi[1]=col2;
    vi[2]=col3;
    read(str,vi);
  }
  void read(const std::string &str, int col1, int col2){
    std::vector<int> vi(2);
    vi[0]=col1; 
    vi[1]=col2;
    read(str,vi);
  }
  void read(const std::string &str, int col){
    std::vector<int> vi(1,col);
    read(str,vi);
  }
  void read(const std::string &str){
    int ncol=colums_in_first_row(str);
    std::vector<int> indices;
    for(int i=0; i<ncol; i++)indices.push_back(i);
    read(str, indices);
  }
  ReadTable(const std::string &str, std::vector<int> indices){
    read(str, indices);
  }
  ReadTable(const std::string &str, int col1, int col2, int col3){
    read(str, col1, col2, col3);
  }
  ReadTable(const std::string &str, int col1, int col2){
    read(str, col1, col2);
  }
  ReadTable(const std::string &str, int col){
    read(str, col);
  }
  ReadTable(const std::string &str){
    read(str);
  }
};

#endif
