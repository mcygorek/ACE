#ifndef AUGMENTED_DENSITY_MATRIX_MPS_DEFINED_H
#define AUGMENTED_DENSITY_MATRIX_MPS_DEFINED_H

#include "Propagator.hpp"
#include <Eigen/Dense>
#include "MPS.hpp"
#include "InfluenceFunctional_Vector.hpp"
#include "RankCompressor.hpp"
//#include "Tensor.hpp"

namespace ACE{

class AugmentedDensityMatrix_MPS{
public:
  Eigen::MatrixXcd rho;
  int n_max;
  MPS ten;

  inline int get_NL()const{return ten.get_dim(0);}
  inline int get_n_max()const{return n_max;}
  inline int get_Ngrps2()const{
    if(ten.get_rank()<2)return ten.get_dim(0);
    else return ten.get_dim(1);
  }

  void check_dimensions()const;

  virtual void print_status(std::ostream &os=std::cout)const;

  void update_rho(int step);

  void propagate(Propagator &prop, const InfluenceFunctional_Vector &IF, 
                 double t, double dt, int step, RankCompressor *compressor,
                 bool use_symmetric_Trotter);
  
  AugmentedDensityMatrix_MPS(int n_max_, int Ngrps, const Eigen::MatrixXcd &rho_);
};

}
#endif
