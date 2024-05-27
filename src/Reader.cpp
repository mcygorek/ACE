#include "PCH.hpp"
#include "ReadExpression.hpp"
#include "ReaderBasics.hpp"
#include "Reader.hpp"

namespace ACE{

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

double readDouble(const std::string &str, const std::string &field){
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

size_t readSizeT(const std::string &str, const std::string &field, size_t max){
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

std::string Matrix_as_parameter(const Eigen::MatrixXcd &M, double epsilon){
  int d=M.rows();
  std::stringstream ss;
  bool is_first=true;
  ss<<"{";
  for(int i=0; i<d; i++){
    for(int j=0; j<d; j++){
      if(abs(M(i,j))<epsilon)continue;

      double re=M(i,j).real();
      double im=M(i,j).imag();

      if(is_first){
        is_first=false;
      }else{
        ss<<"+";
      }

      ss<<"("<<re;
      if(im>epsilon){
        ss<<"+i*("<<im<<")";
      }
      ss<<")*|"<<i<<"><"<<j<<"|_"<<int_to_string(d);
    }
  }
  if(is_first){ss<<"0*Id_"<<int_to_string(d);}
  ss<<"}";
  return ss.str();
}
Eigen::MatrixXcd round_Matrix(const Eigen::MatrixXcd &M, double epsilon){
  if(epsilon==0.)return M;
  Eigen::MatrixXcd tmp(M.rows(), M.cols());
  for(int r=0; r<tmp.rows(); r++){
    for(int c=0; c<tmp.cols(); c++){
      std::complex<double> e=M(r,c);
      if(abs(e.real())<epsilon)e.real(0.);
      if(abs(e.imag())<epsilon)e.imag(0.);
      tmp(r,c)=e;
    }
  }
  return tmp;
}
}//namespace

