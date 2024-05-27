#include "PCH.hpp"
#include "BinaryReader.hpp"
#include "Reader.hpp"
#include "ReadTable.hpp"
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include "DummyException.hpp"

namespace ACE{

binary_write_int_type const binary_write_int=&binary_write<int>;

int binary_read_int(std::istream &ifs, const std::string context){
  int i;
  ifs.read((char*)&i, sizeof(int));
  if(!ifs.good()){
    std::cerr<<"binary_read_int: cannot read int!";
    if(context!="")std::cerr<<" in context: "<<context;
    std::cerr<<std::endl;
    throw DummyException();
  }
  return i;
}

void binary_write_fixedSizeString(std::ostream &ofs, int sz, const std::string &str){
  if(sz!=str.size()){
    std::cerr<<"binary_write_fixedSizeString: sz!=str.size() ("<<sz<<" vs. '"<<str<<"'!"<<std::endl;
    throw DummyException();
  }
  ofs.write((char*)str.c_str(), sizeof(char)*sz);
}
std::string binary_read_fixedSizeString(std::istream &ifs, int sz, const std::string context){
  if(sz<1){
    std::cerr<<"binary_read_fixedSizeString: sz<1";
    if(context!="")std::cerr<<" in context: "<<context;
    std::cerr<<std::endl;
    throw DummyException();
  }
  char *cs=new char[sz+1];
  cs[sz]=0;
  ifs.read((char*)cs, sizeof(char)*sz);
  if(!ifs.good()){
    std::cerr<<"binary_read_fixedSizeString: cannot read fixedSizeString!";
    if(context!="")std::cerr<<" in context: "<<context;
    std::cerr<<std::endl;
    throw DummyException();
  }
  std::string str(cs);
  delete[] cs;
  return str;
}

void binary_write_string(std::ostream &ofs, const std::string &str){
  int sz=str.size();
  binary_write_int(ofs, sz); 
  binary_write_fixedSizeString(ofs, sz, str);
}
std::string binary_read_string(std::istream &ifs, const std::string context){
  int sz=binary_read_int(ifs, context);
  if(sz<0||sz>1024){
    std::cerr<<"binary_read_string: string larger than 1024 characters";
    if(context!="")std::cerr<<" in context '"<<context<<"'";
    std::cerr<<"!"<<std::endl;
    throw DummyException();
  }
  return binary_read_fixedSizeString(ifs, sz, context);
}


void binary_write_EigenMatrixXd(std::ostream &ofs, const Eigen::MatrixXd &M){
  binary_write(ofs, (int)M.rows());
  binary_write(ofs, (int)M.cols());
  double *d=new double[M.rows()*M.cols()];
  for(int r=0; r<M.rows(); r++){
    for(int c=0; c<M.cols(); c++){
      d[r*M.cols()+c]=M(r,c);
    }
  }
  ofs.write((char*)d, sizeof(double)*M.rows()*M.cols());
  delete[] d;
}
Eigen::MatrixXd binary_read_EigenMatrixXd(std::istream &ifs, const std::string context){
  int rows=binary_read_int(ifs, context==""?"":context+" (rows)");
  int cols=binary_read_int(ifs, context==""?"":context+" (cols)");

  Eigen::MatrixXd M(rows, cols);
  double *d=new double[rows*cols];
  ifs.read((char*)d, sizeof(double)*rows*cols);
  for(int r=0; r<M.rows(); r++){
    for(int c=0; c<M.cols(); c++){
      M(r,c)=d[r*M.cols()+c];
    }
  }
  delete[] d;

  if(!ifs.good()){
    std::cerr<<"binary_read_EigenMatrixXd: cannot read EigenMatrixXd!";
    if(context!="")std::cerr<<" in context: "<<context;
    std::cerr<<std::endl;
    throw DummyException();
  }
  return M;
}
void binary_write_EigenMatrixXcd(std::ostream &ofs, const Eigen::MatrixXcd &M){
  binary_write(ofs, (int)M.rows());
  binary_write(ofs, (int)M.cols());
  std::complex<double> *d=new std::complex<double>[M.rows()*M.cols()];
  for(int r=0; r<M.rows(); r++){
    for(int c=0; c<M.cols(); c++){
      d[r*M.cols()+c]=M(r,c);
    }
  }
  ofs.write((char*)d, sizeof(std::complex<double>)*M.rows()*M.cols());
  delete[] d;
}
Eigen::MatrixXcd binary_read_EigenMatrixXcd(std::istream &ifs, const std::string context){
  int rows=binary_read_int(ifs, context==""?"":context+" (rows)");
  int cols=binary_read_int(ifs, context==""?"":context+" (cols)");

  Eigen::MatrixXcd M(rows, cols);
  std::complex<double> *d=new std::complex<double>[rows*cols];
  ifs.read((char*)d, sizeof(std::complex<double>)*rows*cols);
  for(int r=0; r<M.rows(); r++){
    for(int c=0; c<M.cols(); c++){
      M(r,c)=d[r*M.cols()+c];
    }
  }
  delete[] d;

  if(!ifs.good()){
    std::cerr<<"binary_read_EigenMatrixXcd: cannot read EigenMatrixXcd!";
    if(context!="")std::cerr<<" in context: "<<context;
    std::cerr<<std::endl;
    throw DummyException();
  }
  return M;
}

std::string read_first_bytes(const std::string & fname, int n){
  std::string str="";
  std::ifstream ifs(fname.c_str());
  if(!ifs.good())return str;

  for(int i=0; i<n; i++){
    char c=ifs.get();
//std::cout<<"i: "<<i<<" c="<<c<<" str='"<<str<<"'"<<std::endl;
    if(!ifs.good())return str;
    str+=c;
  }
  return str;
}

void write_MatrixXcd(const Eigen::MatrixXcd &A, const std::string &filename){
  std::ofstream ofs(filename.c_str());
  ofs<<"#MatrixXcd: "<<A.rows()<<"x"<<A.cols()<<std::endl;
  for(int r=0; r<A.rows(); r++){
    for(int c=0; c<A.cols(); c++){
      ofs<<A(r,c).real()<<" "<<A(r,c).imag();
      if(c<A.cols()-1){
        ofs<<" ";
      }else{
        ofs<<std::endl;
      }
    }
  }
}
Eigen::MatrixXcd read_MatrixXcd(const std::string &filename){
  ReadTable tab(filename);
  if(tab.table.size()<1){
    std::cerr<<"read_MatrixXcd '"<<filename<<"': fewer than 1 row!"<<std::endl; 
    throw DummyException();
  }
  if(tab.table[0].size()<2){
    std::cerr<<"read_MatrixXcd '"<<filename<<"': fewer than 2 columns!"<<std::endl; 
    throw DummyException();
  }
  if(tab.table[0].size()%2!=0){
    std::cerr<<"read_MatrixXcd '"<<filename<<"': columns not divisible by 2!"<<std::endl; 
    throw DummyException();
  }
  Eigen::MatrixXcd A((int)tab.table.size(), (int)(tab.table[0].size()/2));
  for(int r=0; r<(int)tab.table.size(); r++){
    for(int c=0; c<(int)tab.table[0].size()/2; c++){
      A(r,c)=std::complex<double>(tab[r][2*c], tab[r][2*c+1]);
    }
  }
  return A;
}

}
