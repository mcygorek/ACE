#ifndef READ_TABLE_DEFINED_H
#define READ_TABLE_DEFINED_H

#include <vector>
#include <iosfwd>
//#include "Reader.hpp"

namespace ACE{

class ReadTable{
public:
  std::vector<std::vector<double> > table;

  inline size_t size()const{
    return table.size();
  }

  std::vector<double> & operator[] (int i);
  
  inline operator std::vector<std::vector<double> > () const{
    return table;
  }
  static int colums_in_first_row(const std::string &fname);
  
  void read(const std::string &str, std::vector<int> indices);
  void read(const std::string &str, int col1, int col2, int col3);
  void read(const std::string &str, int col1, int col2);
  void read(const std::string &str, int col);
  void read(const std::string &str);
  
  inline ReadTable(const std::string &str, std::vector<int> indices){
    read(str, indices);
  }
  inline ReadTable(const std::string &str, int col1, int col2, int col3){
    read(str, col1, col2, col3);
  }
  inline ReadTable(const std::string &str, int col1, int col2){
    read(str, col1, col2);
  }
  inline ReadTable(const std::string &str, int col){
    read(str, col);
  }
  inline ReadTable(const std::string &str){
    read(str);
  }
};

}//namespace
#endif
