#ifndef ACE_ENVIRONMENT_OPERATORS_DEFINED_H
#define ACE_ENVIRONMENT_OPERATORS_DEFINED_H

#include "PassOn.hpp"
#include "SelectIndices.hpp"
#include <iostream>
#include <vector>

namespace ACE{

class EnvironmentOperators{
public:
  std::vector<Eigen::VectorXcd> ops;
  int use_fermion;
 
  inline void swap(EnvironmentOperators &other){
    ops.swap(other.ops);
    std::swap(use_fermion, other.use_fermion);
  }
  inline void clear(){ 
    ops.clear();
    use_fermion=0;
  }
  inline void remove_all_but_first(){ 
    ops.resize(1);
    use_fermion=0;
  }
  inline bool is_used()const{
    return (ops.size()>0);
  }

  inline operator const std::vector<Eigen::VectorXcd> & () const{
    return ops;
  }
  inline size_t size()const{return ops.size();}
  const Eigen::VectorXcd & operator[] (int i)const;

  void set_from_matrices(const std::vector<Eigen::MatrixXcd> & Mvec);
  void join(const EnvironmentOperators &other);
  void join_select_indices(const EnvironmentOperators &other,
                           const SelectIndices &k_right);

  void print_debug(std::ostream &ofs=std::cout)const;

  void set_ill_defined();
  void process_forward(const PassOn & pass_on);
  void process_backward(const PassOn & pass_on);

  void read_binary(std::istream &is);
  void write_binary(std::ostream &os)const;

  EnvironmentOperators() : use_fermion(0){
  }
  EnvironmentOperators(const std::vector<Eigen::MatrixXcd> & Mvec, int fermion=0){
    set_from_matrices(Mvec);
    use_fermion=fermion;
  }
};


}//namespace
#endif
