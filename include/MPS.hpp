#ifndef MPS_DEFINED_H_
#define MPS_DEFINED_H_

#include "Tensor.hpp"
#include "Eigen_fwd.hpp"
#include <iosfwd>
#include <complex>
#include "MPS_Matrix.hpp"

/** Goal:

We want to be able to store (and propagate) a tensor
in the form of a matrix-product state:

A_{i_0, i_1, ... i_{n-1}}   = 
\sum_{d_1,d_2,.. d_{n-1}} a_{i_0}^{0 d_1} a_{i_1}^{d_1d_2} .. a_{i_{n-1}}^{d_{n-1} 0 } 

The decomposition will be performed using sweeps of SVDs.


Every "a" is an object with three different indices

*/

namespace ACE{

template <typename T> class MPS_ScalarType: public Tensor_ScalarType<T>{
public:
  std::vector<MPS_Matrix_ScalarType<T> > a;
 
  inline virtual int get_rank() const{return a.size();}
  inline virtual int get_dim(int i) const{return a[i].dim_i;}
  inline size_t size()const{return a.size();}
  inline MPS_Matrix_ScalarType<T> & back(){return a.back();}
  inline const MPS_Matrix_ScalarType<T> & back() const{return a.back();}

  virtual void resize(const std::vector<int> &list);
  
  void resize(int rank, int dim);
  
  void resize_fill_one(int rank, int dim);
  
  void reduce_last(const Eigen::Matrix<T, Eigen::Dynamic, 1> &v);
  void reduce_first(const Eigen::Matrix<T, Eigen::Dynamic, 1> &v);

  void print_dims(std::ostream &os=std::cout)const;
  void print_dims(const std::string &filename)const;
  
  int get_max_dim()const;
  int get_max_dim_pos()const;
  
  void print_max_dim(std::ostream &os=std::cout)const;
  void print_inner_dims(std::ostream &os=std::cout)const;

  double max_element_abs()const;

  virtual void check_consistency(const std::vector<int> &i_list=std::vector<int>())const;
  
  T operator()(const std::vector<int> &i_list)const;
  inline T operator()(const Tensor_Index &ind)const{
    return operator()(ind.list);
  }

  void multiply_front(const MPS_ScalarType<T> & other);
  void low_end_multiply(const MPS_ScalarType<T> & other);
  void insert(int pos, const MPS_Matrix_ScalarType<T> &b, int n_insert=1);
  
  void write_binary(std::ostream &ofs)const;
  void write_binary(const std::string &filename)const;
  
  virtual void read_binary(std::istream &ifs);
  virtual void read_binary(const std::string &filename);
  
  void copy(const MPS_ScalarType<T> &other);
  
  inline void swap(MPS_ScalarType<T> &other){
    a.swap(other.a);
  }
  inline MPS_ScalarType &operator=(const MPS_ScalarType<T> &other){
    copy(other);
    return *this;
  }
  inline MPS_ScalarType(const MPS_ScalarType<T> &other){
    copy(other);
  }
  inline MPS_ScalarType(){}
  virtual ~MPS_ScalarType(){}

  static MPS_ScalarType<T> get_trivial(int nmax, int dim_i);
};

typedef MPS_ScalarType<std::complex<double> > MPS;
typedef MPS_ScalarType<double> MPS_real;

}//namespace
#endif
