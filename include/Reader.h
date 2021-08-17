#ifndef READER_DEFINED_H
#define READER_DEFINED_H

#include "ReadExpression.h"

namespace Reader{

///returns true is string can be interpreted as double
///may also be more complicated expression within { } 
bool canReadDouble(const std::string &str, double &d){
  if(canReadDouble_strict(str, d)){
    return true;
  }
  if(str.size()>0 && str[0]=='{'){
    Eigen::MatrixXcd M=ReadExpression(str);
    if(M.rows()!=1)return false;
    d=M(0,0).real(); 
    return true;
  }
  return false;
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

int number_of_doubles(const std::string &line){
  std::vector<std::string> toks=tokenize(line);
  for(size_t i=0; i<toks.size(); i++){
    if(!isDouble(toks[i]))return i;
  }
  return toks.size();
}

}//namespace

std::string add_prefix(const std::string & prefix, const std::string &str){
  if(prefix=="")return str;
  else return std::string(prefix+"_"+str);
}

#endif
