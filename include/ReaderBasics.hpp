#pragma once
#ifndef READER_BASICS_DEFINED_H
#define READER_BASICS_DEFINED_H

#include <iosfwd>
#include <string>
//#include <iostream>
#include <fstream>
//#include <sstream>
//#include <string>
#include <vector>
#include <utility>
#include <memory>


/** 
Hard core of "Reader" functions to deal with basic interpretation of 
strings. The actual Reader.hpp uses functions from ReadExpression.hpp, which
in turn accesses more basic functions. These are define here. 
*/


namespace ACE{


///returns true is string can be interpreted as double
bool canReadDouble_strict(const std::string &str, double &d);

inline bool isDouble_strict(const std::string &str){
  double d;
  return canReadDouble_strict(str,d);
}

std::string int_to_string(int d);

std::string double_to_string(double d);

std::string trim(const std::string &str);


//Note: only moves stream forward if string 'str' is found
bool is_next_in_stream(std::stringstream & ss, const std::string & str);

inline bool is_next_in_stream(std::stringstream & ss, char c){
  return is_next_in_stream(ss, std::string(1,c));
}

bool find_matching_brace_block(std::stringstream &ss, std::string &result, char cbegin='(', char cend=')');


std::string find_matching_brace_block_complain(std::stringstream &ss, char cbegin='(', char cend=')');


bool getRelevantLine(std::istream &is, std::string &line);

int getNrRelevantLines(std::istream &is);

inline int getNrRelevantLines(const std::string &fname){
  std::ifstream ifs(fname.c_str());
  return getNrRelevantLines(ifs);
}

void tokenize(const std::string &line, std::vector<std::string> & toks);

inline std::vector<std::string> tokenize(const std::string &line){
  std::vector<std::string> toks;
  tokenize(line, toks);
  return toks;
}

bool getRelevantLineTokens(std::istream &is, std::vector<std::string> & toks);

std::string stringvec_to_string(const std::vector<std::string> &sv);

bool file_exists(const std::string &str);
void print_file_exists(const std::string &str);
void check_file_exists(const std::string &str, bool verbose=false);

std::unique_ptr<std::ifstream> open_file_check(const std::string &str);

}//namespace

#endif
