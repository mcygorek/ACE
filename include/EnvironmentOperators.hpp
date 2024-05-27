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
 
  inline void swap(EnvironmentOperators &other){
    ops.swap(other.ops);
  }
  inline void clear(){ ops.clear(); }
  inline void remove_all_but_first(){ ops.resize(1);}
  inline bool is_used()const{
    return (ops.size()>0);
  }

  inline operator const std::vector<Eigen::VectorXcd> & () const{
    return ops;
  }
  inline size_t size()const{return ops.size();}
  const Eigen::VectorXcd & operator[] (int i)const;

  void set_from_matrices(const std::vector<Eigen::MatrixXcd> & Mvec, int N_mode);
  void join(const EnvironmentOperators &other);
  void join_select_indices(const EnvironmentOperators &other,
                           const SelectIndices &k_right);

  void print_debug(std::ostream &ofs=std::cout)const;

  void set_ill_defined();
  void process_forward(const PassOn & pass_on);
  void process_backward(const PassOn & pass_on);

  void read_binary(std::istream &is);
  void write_binary(std::ostream &os)const;
};


}//namespace
#endif
