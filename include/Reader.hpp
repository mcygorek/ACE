#pragma once
#ifndef READER_DEFINED_H
#define READER_DEFINED_H

#include "ReadExpression.hpp"
#include "ReaderBasics.hpp"

namespace ACE{

///returns true is string can be interpreted as double
///may also be more complicated expression within { } 
bool canReadDouble(const std::string &str, double &d);

bool isDouble(const std::string &str);

double readDouble(const std::string &str, const std::string &field="");

size_t readSizeT(const std::string &str, const std::string &field="", size_t max=0);

int number_of_doubles(const std::string &line);

inline std::string add_prefix(const std::string & prefix, const std::string &str){
  if(prefix=="")return str;
  else return std::string(prefix+"_"+str);
}

std::string Matrix_as_parameter(const Eigen::MatrixXcd &M, double epsilon=1e-10);

Eigen::MatrixXcd round_Matrix(const Eigen::MatrixXcd &M, double epsilon);

}//namespace
#endif
