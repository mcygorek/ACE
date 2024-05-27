#ifndef IF_OD_DICTIONARY_DEFINED_H
#define IF_OD_DICTIONARY_DEFINED_H

#include "MPS_Matrix.hpp"
#include "ReadPT_struct.hpp"
#include "DiagBB.hpp"
#include <vector>

namespace ACE{
template <typename T> class MPS_ScalarType;



/**  maps (alpha, tilde{alpha}) \to beta; store identical only once */
class IF_OD_Dictionary{
public:
  int N;  
  int reduced_dim;
  //Mapping. Note: if -1: contribution is zero
  std::vector<int> beta;


  inline void swap(IF_OD_Dictionary &other){
    int tmp=N; N=other.N; other.N=tmp;
    tmp=reduced_dim; reduced_dim=other.reduced_dim; other.reduced_dim=tmp;
    beta.swap(other.beta);
  }
  inline int get_N()const{ return N; }
  inline int get_NL()const{ return get_N()*get_N(); }
  inline int get_NL2()const{ return get_NL()*get_NL(); }
  inline int get_reduced_dim()const{ return reduced_dim; }

  void calculate_reduced_dim();
  bool is_diagonal()const;
  
  inline int operator()(int i)const{ return beta[i];}

  bool operator==(const IF_OD_Dictionary &other)const;
  inline bool operator!=(const IF_OD_Dictionary &other)const{
    return !operator==(other);
  }

  void print_beta(std::ostream &os=std::cout)const;

  void set_default(int n);
  void set_default_diag(int n);
  void set_trivial(int n);

  std::vector<std::vector<int> > get_reverse_beta()const;
  
  /** 
    For joining two IF_ODs with possibly different coupling structure:
    Find smallest common dictionary. 
    Returns true if it has changed
  */
  bool join(const IF_OD_Dictionary &other);
  
  template <typename T> 
  void detect(const MPS_Matrix_ScalarType<T> &a, double zero=1e-12);
  
  template <typename T>
  void detect(const MPS_ScalarType<T> &m, double zero=1e-12);
  
  double get_keep_weight()const;

  void read_binary(std::istream &ifs);
  
  void write_binary(std::ostream &ofs)const;
  
  //Assume $\mathcal{H_E}\to \mathcal{H_X}\otimes \mathcal{H_E}$
  void expand_space_front(int Nfront);
  void expand_space_back(int Nback);
  void expand(const ReadPT_struct &read_PT, bool verbose=true);

  void expand_DiagBB(const DiagBB &diagBB);

  //reduce from dense to sparse
  template <typename T>
  void reduce_MPS_Matrix(MPS_Matrix_ScalarType<T> &M)const;
  //expande from sparse to dense
  template <typename T>
  void expand_MPS_Matrix(MPS_Matrix_ScalarType<T> &M)const;
  template <typename T>
  //translate from 'other' sparse form to '*this' sparse form
  void translate_MPS_Matrix(MPS_Matrix_ScalarType<T> &M, const IF_OD_Dictionary &other)const;
 
  //undo group decomposition for diagonal coupling. Modify both dict and M.
  template <typename T>
  void expand_all_diagonal(MPS_Matrix_ScalarType<T> &M);
 
  void copy(const IF_OD_Dictionary &other);
  
  inline IF_OD_Dictionary & operator=(const IF_OD_Dictionary &other){
    copy(other);
    return *this;
  }
  inline IF_OD_Dictionary (const IF_OD_Dictionary &other){
    copy(other);
  }
  template <typename T>
  IF_OD_Dictionary(const MPS_Matrix_ScalarType<T> &a, double zero=1e-12){
    detect(a, zero);
  }
  template <typename T>
  IF_OD_Dictionary(const MPS_ScalarType<T> &m, double zero=1e-12){
    detect(m, zero);
  }
  inline IF_OD_Dictionary(int def=2){
    set_default(def);
  }
};

}//namespace
#endif
