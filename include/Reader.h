#ifndef READER_DEFINED_H
#define READER_DEFINED_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>

namespace Reader{


///returns true is string can be interpreted as double
bool canReadDouble(const std::string &str, double &d){
  std::stringstream ss(str);
  ss>>d;
  if(ss.fail()) return false;
  else return true; 
}
bool isDouble(const std::string &str){
  double d;
  return canReadDouble(str,d);
}
double readDouble(const std::string &str, const std::string &field=""){
  double d;
  if(!canReadDouble(str, d)){
    if(field==""){
      std::cerr<<"Cannot convert '"<<str<<"' into double!"<<std::endl;
      exit(1);
    }else{
      std::cerr<<"Cannot convert '"<<str<<"' into double in context '"<<field<<"'!"<<std::endl;
      exit(1);
    }
  }
  return d;
}
size_t readSizeT(const std::string &str, const std::string &field="", size_t max=0){
  int i=round(readDouble(str,field));
  if(i<0){std::cerr<<"'"<<field<<"': '"<<str<<"' <0!"<<std::endl; exit(1);}
  if(max>0 && (size_t)i>=max){
    std::cerr<<"'"<<field<<"': '"<<str<<"' should be smaller than "<<max<<"!"<<std::endl; exit(1);}
  return i;
}

std::string trim(const std::string &str){
  int nr_begin=0;
  int nr_end=0;
  char c;
  for(nr_begin=0; nr_begin<str.length(); nr_begin++){
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
bool is_next_in_stream(std::stringstream & ss, char c){
  return is_next_in_stream(ss, std::string(1,c));
}
bool find_matching_brace_block(std::stringstream &ss, std::string &result, char cbegin='(', char cend=')'){

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
std::string find_matching_brace_block_complain(std::stringstream &ss, char cbegin='(', char cend=')'){
  std::string result;
  if(!find_matching_brace_block(ss, result, cbegin, cend)){
    std::cerr<<"Can't find matching brace block with delimiters '"<<cbegin<<"', '"<<cend<<"'!"<<std::endl;
    exit(1);
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
std::vector<std::string> tokenize(const std::string &line){
  std::vector<std::string> toks;
  tokenize(line, toks);
  return toks;
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

}

#endif
