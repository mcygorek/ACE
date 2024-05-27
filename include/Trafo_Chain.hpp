#ifndef ACE_TRAFO_CHAIN_DEFINED_H
#define ACE_TRAFO_CHAIN_DEFINED_H

#include <vector>
#include <Eigen/Dense>
#include "MPS.hpp"
#include "ModePropagatorGenerator.hpp"
#include "Tensor.hpp"
#include "BinaryReader.hpp"

namespace ACE{

Eigen::MatrixXd pseudoinverse(const Eigen::MatrixXd &M);

class Trafo_Chain{
public:
  std::vector<Eigen::MatrixXd> T;
  std::vector<Eigen::MatrixXd> Tinv;

  inline size_t size()const{return T.size();}
  inline int lastdim()const{
    if(T.size()<1)return 0;
    return T.back().rows();
  }

  void print_info(std::ostream &ofs=std::cout);

  void add_low_to_high(const Eigen::MatrixXd &R);

  void add_high_to_low(const Eigen::MatrixXd &L);

  std::vector<int> get_dims()const;

  void check_compatible(const Trafo_Chain &other)const;
  
  Eigen::MatrixXd overlap_matrix(const Trafo_Chain &other, bool skiplastmultiply=false)const;
  
  double overlap(const Trafo_Chain &other)const;

  void combine(const Trafo_Chain &other, int max_add=-1);

  void orthogonalize_single(int i, int j);

  void orthogonalize();

  //Perform SVD on all but the last T:
  void compress_all_but_last_T(double epsilon);
  
  //Perform SVD on all but the last Tinv:
  void compress_all_but_last_Tinv(double epsilon); 
  
  void compress_weight_to_Tinv(double epsilon);  //make orthogonal. Residual SVD to last Tinv.

  void compress_weight_sym(double epsilon);

  void SVD_sweep_T_to_Tinv(double epsilon);
  
  //Perform SVD on all but the last Tinv:
  void compress(double epsilon);

  void SVD_orthogonalize(double epsilon);

  void Eigen_orthogonalize(double epsilon);

  void add_ortho(Trafo_Chain other, double epsilon, double epsilon2=0, int max_add=-1);

  void add_ortho2(Trafo_Chain other, double epsilon, int max_add=-1);

  void combine_Tinv_T(double epsilon);

  void read(const std::string &fname);
  
  void write(const std::string &fname);
  
  inline Trafo_Chain(const std::string &fname){
    read(fname);
  }
  inline Trafo_Chain(){}
  inline virtual ~Trafo_Chain(){}
};
}//namespace
#endif
