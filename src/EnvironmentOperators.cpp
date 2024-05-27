#include "EnvironmentOperators.hpp"
#include "LiouvilleTools.hpp"
#include "otimes.hpp"
#include <iostream>
#include "BinaryReader.hpp"
#include "DummyException.hpp"

namespace ACE{

const Eigen::VectorXcd & EnvironmentOperators::operator[] (int i)const{
  if(i<0||i>=ops.size()){
    std::cerr<<"EnvironmentOperators::operator[]: out of bounds: i<0||i>=ops.size() ("<<i<<" vs. "<<ops.size()<<")!"<<std::endl;
    throw DummyException();
  }
  return ops[i];
}

void EnvironmentOperators::set_from_matrices(const std::vector<Eigen::MatrixXcd> & Mvec, int N_mode){
  /*
  if(Mvec.size()<1){
    ops.clear();
    return;
  }
  */
  ops.resize(Mvec.size()+1);
  ops[0]=H_Matrix_to_L_Vector(Eigen::MatrixXcd::Identity(N_mode,N_mode));
  for(int i=0; i<(int)Mvec.size(); i++){
    if(Mvec[i].rows()!=N_mode || Mvec[i].cols()!=N_mode){
      std::cerr<<"EnvironmentOperators::set_from_matrices: Mvec[i].rows()!=N || Mvec[i].cols()!=N!"<<std::endl;
      throw DummyException();
    }
    ops[i+1]=H_Matrix_to_L_Vector(Mvec[i].transpose());
  }
}

void EnvironmentOperators::join(const EnvironmentOperators &other){
//std::cout<<"EnvOps::join: size()="<<size()<<" other.size()="<<other.size()<<std::endl;
  if(size()<1 || other.size()<1){
    ops.clear();
    return;
  }
  if(other.size()==1){
    for(size_t i=0; i<size(); i++){
      ops[i]=Vector_otimes(ops[i], other[0]);
    }
    return;
  }

  if(size()==1){
    std::vector<Eigen::VectorXcd> new_ops(other.size());
    for(size_t i=0; i<other.size(); i++){
      new_ops[i]=Vector_otimes(ops[0], other[i]);
    }
    ops.swap(new_ops);
    return;
  }
 
  if(size()!=other.size()){
    std::cerr<<"EnvironmentOperators::join: size()!=other.size() ("<<size()<<" vs. "<<other.size()<<")!"<<std::endl;
    throw DummyException();
  }

  std::vector<Eigen::VectorXcd> new_ops(size()); 
  new_ops[0]=Vector_otimes(ops[0], other[0]);
  for(size_t i=1; i<size(); i++){
    new_ops[i] = Vector_otimes(ops[i], other[0])  
               + Vector_otimes(ops[0], other[i]);
  }
  ops.swap(new_ops);
}

void EnvironmentOperators::join_select_indices(
        const EnvironmentOperators &other, const SelectIndices &k_right){

//std::cout<<"EnvOps::join_select_indices: size()="<<size()<<" other.size()="<<other.size()<<std::endl;
  if(size()<1 || other.size()<1){
    ops.clear();
    return;
  }
  if(other.size()==1){
    for(size_t i=0; i<size(); i++){
      ops[i]=k_right.Vector_otimes(ops[i], other[0]);
    }
    return;
  }

  if(size()==1){
    std::vector<Eigen::VectorXcd> new_ops(other.size());
    for(size_t i=0; i<other.size(); i++){
      new_ops[i]=k_right.Vector_otimes(ops[0], other[i]);
    }
    ops.swap(new_ops);
    return;
  }
 
  if(size()!=other.size()){
    std::cerr<<"EnvironmentOperators::join: size()!=other.size() ("<<size()<<" vs. "<<other.size()<<")!"<<std::endl;
    throw DummyException();
  }

  std::vector<Eigen::VectorXcd> new_ops(size()); 
  new_ops[0]=k_right.Vector_otimes(ops[0], other[0]);
  for(size_t i=1; i<size(); i++){
    new_ops[i] = k_right.Vector_otimes(ops[i], other[0])  
               + k_right.Vector_otimes(ops[0], other[i]);
  }
  ops.swap(new_ops);
}

void EnvironmentOperators::print_debug(std::ostream &ofs)const{
  ofs<<"size: "<<ops.size()<<" [";
  for(size_t i=0; i<ops.size(); i++)ofs<<" "<<ops[i].rows();
  ofs<<" ]";
}

void EnvironmentOperators::set_ill_defined(){
  if(ops.size()>0){
    ops[0]=Eigen::VectorXcd::Ones(1);
  }
  for(int i=1; i<(int)ops.size(); i++){
    ops[i]=Eigen::VectorXcd::Ones(1)*(1./0.);
  }
}

void EnvironmentOperators::process_forward(const PassOn & pass_on){
  for(size_t o=0; o<ops.size(); o++){
    Eigen::VectorXcd & op=ops[o];
    if(pass_on.P.cols()!=op.rows()){
      std::cerr<<"EnvironmentOperators::process_forward: pass_on.P.cols()="<<pass_on.P.cols()<<"!=ops["<<o<<"].rows()="<<op.rows()<<"!"<<std::endl;
      throw DummyException();
    }
    op = pass_on.P * op;
  }
}
void EnvironmentOperators::process_backward(const PassOn & pass_on){
  for(Eigen::VectorXcd & op : ops){
    if(op.rows()!=pass_on.Pinv.cols()){
      std::cerr<<"EnvironmentOperators::process_backward: op.rows()!=pass_on.Pinv.cols() ("<<op.rows()<<" vs. "<<pass_on.Pinv.cols()<<std::endl;
      throw DummyException();
    }
    op = pass_on.Pinv * op;
  }
}

void EnvironmentOperators::read_binary(std::istream &is){
  int sz=binary_read_int(is, "EnvOps");
  if(sz<0){
    std::cerr<<"EnvironmentOperators::read_binary: sz<0!"<<std::endl;
    throw DummyException();
  } 
  ops.resize(sz);
  for(int i=0; i<sz; i++){
    ops[i]=binary_read_EigenMatrixXcd(is, "EnvOps");
  }
}
void EnvironmentOperators::write_binary(std::ostream &os)const{
  int sz=ops.size();
  binary_write_int(os, sz);
  for(int i=0; i<sz; i++){
    binary_write_EigenMatrixXcd(os, ops[i]);
  }
} 

}//namespace
