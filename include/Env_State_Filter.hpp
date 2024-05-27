#ifndef ACE_ENV_STATE_FILTER_DEFINED_H
#define ACE_ENV_STATE_FILTER_DEFINED_H

#include "MPS.hpp"
#include "Eigen_fwd.hpp"

namespace ACE{
class Parameters;

template <typename Scalar_Type>
Eigen::Matrix<Scalar_Type, Eigen::Dynamic, Eigen::Dynamic>  order_trafo(
                         Eigen::Matrix<Scalar_Type, Eigen::Dynamic, 1> &q);


class Env_State_Filter{
public:
  bool mean_field;
  bool no_rotate_first; 
  int nr_MF;

  std::vector<std::pair<int,int> > selected;
  int dim1, dim2;
  std::vector<std::pair<int,int> > last_selected;
  int last_dim1, last_dim2;

  void rotate_identity_to_first(std::vector<MPS_Matrix_real> &a, 
                 std::vector<std::vector<Eigen::VectorXcd> > &env_ops, int n);

  void preprocess(std::vector<MPS_Matrix_real> &a1,
                  std::vector<std::vector<Eigen::VectorXcd> > &env_ops1, 
                  int n1,
                  std::vector<MPS_Matrix_real> &a2,
                  std::vector<std::vector<Eigen::VectorXcd> > &env_ops2, 
                  int n2);

  void filter(std::vector<MPS_Matrix_real> &a, 
             std::vector<std::vector<Eigen::VectorXcd> > &env_ops, 
             int n);

  void setup(Parameters &param);
  
  Env_State_Filter(Parameters &param);
  Env_State_Filter();
  inline virtual ~Env_State_Filter(){}
};

}//namespace
#endif
