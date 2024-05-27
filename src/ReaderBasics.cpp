#include "PCH.hpp"
#include "ReaderBasics.hpp"
#include "DummyException.hpp"

#include <iostream>
#include <sstream>
#include <memory>

/** 
Hard core of "Reader" functions to deal with basic interpretation of 
strings. The actual Reader.h uses functions from ReadExpression.h, which
in turn accesses more basic functions. These are define here. 
*/

namespace ACE{


///returns true is string can be interpreted as double
bool canReadDouble_strict(const std::string &str, double &d){
  if(str=="inf"){d=1./0.; return true;}
  std::stringstream ss(str);
  ss>>d;
  if(ss.fail()) return false;
  else return true; 
}

std::string int_to_string(int d){
  std::stringstream ss; ss<<d;
  return ss.str();
}
std::string double_to_string(double d){
  std::stringstream ss; ss<<d;
  return ss.str();
}

std::string trim(const std::string &str){
  int nr_begin=0;
  int nr_end=0;
  for(nr_begin=0; nr_begin<(int)str.length(); nr_begin++){
    if(!isspace(str[nr_begin]))break;
  }
  if(nr_begin==(int)str.length())return ""; 
  for(nr_end=str.length(); nr_end>0; nr_end--){
    if(!isspace(str[nr_end-1]))break;
  }
  return str.substr(nr_begin, nr_end-nr_begin);
}


//Note: only moves stream forward if string 'str' is found
bool is_next_in_stream(std::stringstream & ss, const std::string & str){
  if(!ss.good()){return false;}
  int tell=ss.tellg();
  for(int i=0; i<(int)str.length(); i++){
    char c=ss.get();
    if( (!ss.eof()) && isspace(c) ){i--; continue;}
    if(ss.eof() || c!=str[i]){
      ss.clear(); 
      ss.seekg(tell, std::ios_base::beg); 
      return false;
    }
  }
  return true; 
}

bool find_matching_brace_block(std::stringstream &ss, std::string &result, char cbegin, char cend){

  std::stringstream ss2;
  bool foundfirst=false;
  int nr_open=0;
  char c;
  while(c=ss.get(), !ss.eof()){
    if(!foundfirst){
      if(isspace(c))continue;
      if(c==cbegin){foundfirst=true; continue;}
      return false;
    }
    if(c==cbegin)nr_open++;
    if(c==cend){
      if(nr_open==0){
        result=ss2.str();
        return true;
      }else nr_open--;
    }
    ss2<<c;
  }
  return false;
}
std::string find_matching_brace_block_complain(std::stringstream &ss, char cbegin, char cend){
  std::string result;
  if(!find_matching_brace_block(ss, result, cbegin, cend)){
    std::cerr<<"Can't find matching brace block with delimiters '"<<cbegin<<"', '"<<cend<<"'!"<<std::endl;
    throw DummyException();
  }
  return result;
}

bool getRelevantLine(std::istream &is, std::string &line){
  while(getline(is,line)){
    for(size_t i=0; i<line.length(); i++){
      if(line[i]==' '||line[i]=='\t'||line[i]=='\n')continue;
      if(line[i]=='#')break;
      return true;
    }
  }
  return false;
}
int getNrRelevantLines(std::istream &is){
  int counter=0;
  std::string str;
  while(getRelevantLine(is, str))counter++;
  return counter;
}

void tokenize(const std::string &line, std::vector<std::string> & toks){
  toks.clear();
  std::string str="";
  for(size_t i=0; i<line.length(); i++){
    if(line[i]=='#'){
      if(str!="")toks.push_back(str);
      return;
    }else if(line[i]==' '||line[i]=='\t'||line[i]=='\n'){
      if(str=="")continue;
      toks.push_back(str);
      str="";
    }else{
      str.push_back(line[i]);
    }    
  }
  if(str!="")toks.push_back(str);
}

bool getRelevantLineTokens(std::istream &is, std::vector<std::string> & toks){
  std::string line;
  if(!getRelevantLine(is,line))return false;
  tokenize(line, toks);
  return true;
}
std::string stringvec_to_string(const std::vector<std::string> &sv){
  if(sv.size()<1)return "";
  std::stringstream ss;
  ss<<sv[0];
  for(size_t i=1; i<sv.size(); i++)ss<<" "<<sv[i];
  return ss.str();
}

bool file_exists(const std::string &str){
  std::ifstream ifs(str);
  if(ifs.good()){
    return true;
  }else{
    return false;
  }
}
void print_file_exists(const std::string &str){
  if(file_exists(str)){
    std::cout<<"File '"<<str<<"' exists."<<std::endl;
  }else{
    std::cout<<"File '"<<str<<"' does not exist."<<std::endl;
  }
}
void check_file_exists(const std::string &str, bool verbose){
  if(!file_exists(str)){
    std::cerr<<"Cannot open file '"<<str<<"'!"<<std::endl;
    throw DummyException();
  }
  if(verbose){
    std::cout<<"File '"<<str<<"' exists."<<std::endl;
  }
}
std::unique_ptr<std::ifstream> open_file_check(const std::string &str){
  std::unique_ptr<std::ifstream> ifp(new std::ifstream(str.c_str()));
  if(!ifp|| !(*ifp)){
    std::cerr<<"Cannot open file '"<<str<<"'!"<<std::endl;
    throw DummyException();
  }
  return ifp;
}

}//namespace

