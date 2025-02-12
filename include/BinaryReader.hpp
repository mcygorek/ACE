#pragma once
#ifndef ACE_BINARY_READER_DEFINED_H
#define ACE_BINARY_READER_DEFINED_H

#include "Reader.hpp"
#include <fstream>
#include "DummyException.hpp"

namespace ACE{

template <typename T> void binary_write(std::ostream &ofs, const T &t){
  ofs.write((char*)&t, sizeof(T));
}
template <typename T> T binary_read(std::istream &ifs, const std::string &context=""){
  T t;
  ifs.read((char*)&t, sizeof(T));
  if(!ifs.good()){
    std::cerr<<"binary_read failed!";
    if(context!=""){std::cerr<<"context: \""<<context<<"\"";}
    std::cerr<<std::endl;
    throw DummyException();
  }
  return t;
}
typedef void (*binary_write_int_type)(std::ostream &ofs, const int &i);


extern binary_write_int_type const binary_write_int;
extern int binary_read_int(std::istream &ifs, const std::string context="");

//extern binary_write_vector_int_type const binary_write_vector_int;
//extern int binary_read_vector_int(std::istream &ifs, const std::string context="");


extern void binary_write_fixedSizeString(std::ostream &ofs, int sz, const std::string &str);

extern std::string binary_read_fixedSizeString(std::istream &ifs, int sz, const std::string context="");

extern void binary_write_string(std::ostream &ofs, const std::string &str);

extern std::string binary_read_string(std::istream &ifs, const std::string context="");


extern void binary_write_EigenMatrixXd(std::ostream &ofs, const Eigen::MatrixXd &M);
inline void binary_write_EigenMatrixXd(const std::string &fname, const Eigen::MatrixXd &M){
  std::ofstream ofs(fname);
  binary_write_EigenMatrixXd(ofs, M);
}

extern Eigen::MatrixXd binary_read_EigenMatrixXd(std::istream &ifs, const std::string context="");
inline Eigen::MatrixXd binary_read_EigenMatrixXd(const std::string &fname){
  std::unique_ptr<std::ifstream> ifs = open_file_check(fname);
  return binary_read_EigenMatrixXd(*(ifs.get()), fname);
}

extern void binary_write_EigenMatrixXcd(std::ostream &ofs, const Eigen::MatrixXcd &M);
inline void binary_write_EigenMatrixXcd(const std::string &fname, const Eigen::MatrixXcd &M){
  std::ofstream ofs(fname);
  binary_write_EigenMatrixXcd(ofs, M);
}


extern Eigen::MatrixXcd binary_read_EigenMatrixXcd(std::istream &ifs, const std::string context="");
inline Eigen::MatrixXcd binary_read_EigenMatrixXcd(const std::string &fname){
  std::unique_ptr<std::ifstream> ifs = open_file_check(fname);
  return binary_read_EigenMatrixXcd(*(ifs.get()), fname);
}


extern std::string read_first_bytes(const std::string & fname, int n);


//read and write MatrixXcd as human readable text
void write_MatrixXcd(const Eigen::MatrixXcd &A, const std::string &filename);
Eigen::MatrixXcd read_MatrixXcd(const std::string &filename);


}//namespace
#endif
