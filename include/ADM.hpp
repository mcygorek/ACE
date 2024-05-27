#ifndef AUGMENTED_DENSITY_MATRIX_DEFINED_H
#define AUGMENTED_DENSITY_MATRIX_DEFINED_H

#include "RankCompressor.hpp"
#include "Propagator.hpp"
#include "InfluenceFunctional.hpp"
#include "Tensor_Dense.hpp"

namespace ACE{

class AugmentedDensityMatrix{
public:
  Eigen::MatrixXcd rho;
  Tensor_Dense ten;

  inline int get_NL()const{return ten.get_dim(0);}
  inline int get_n_max()const{return ten.get_rank();}
  inline int get_Ngrps2()const{
    if(ten.get_rank()<2)return ten.get_dim(0);
    else return ten.get_dim(1);
  }

  void check_dimensions()const;
  
  virtual void print_status(std::ostream &os=std::cout)const;

  void update_rho(int step);

  void propagate(Propagator &prop, const InfluenceFunctional &IF, 
            double t, double dt, int step, RankCompressor *compressor,
            bool use_symmetric_Trotter);

  AugmentedDensityMatrix(int n_max, int Ngrps, const Eigen::MatrixXcd &rho_);

};
}
#endif
