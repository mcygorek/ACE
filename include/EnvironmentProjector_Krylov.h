#ifndef ACE_ENVIRONMENT_PROJECTOR_KRYLOV
#define ACE_ENVIRONMENT_PROJECTOR_KRYLOV

#include "EnvironmentProjector.h"

/* To find reduced basis: 
Calculate which env. states (Liouville space)
can be reached with a given trial (e.g. initial) bath states 
by applying the propagator for all possible combinations of 
(\tilde{\alpha},\alpha). 

NOTE: Instead of the propagator, it may be wiser to use its generator, i.e.
the Liouvillian. This will make the analysis independent of the time step.
One the other hand, strongly off-resonant modes will sustain large 
imagniary eigenvalues.

Then, use those as trial states for yet another generation of states.
Orthonormalize wrt. existing states
(Maybe: symmetric orthonormalization wrt. states of same generation).
Store adjacency matrix and use something like PageRank to determine 
importance.
*/

namespace ACE{

class EnvironmentProjector_Krylov: public EnvironmentProjector{
public:

  virtual Eigen::MatrixXcd get_Q(const Eigen::MatrixXcd &Mdis, 
                                 const Eigen::VectorXcd &bath_init) override{

     int ML=get_vector_dim_min(bath_init, 1, "EnvironmentProjector_Schur::get_Q: bath_init");
     int NL=Mdis.rows()/ML;
     check_matrix_square_eq(Mdis, NL*ML, "EnvironmentProjector_Schur::get_Q: Mdis");
     check_at_least(NL, 4, "EnvironmentProjector_Schur::get_Q: NL");
     

     Eigen::MatrixXcd Q=Eigen::MatrixXcd::Identity(ML,ML);
     if(compr.has_effect()){


//       std::cout<<"max_diff_from_ortho(Q): "<<max_diff_from_ortho(Q)<<std::endl;
     }
     return Q;
  }
  virtual Eigen::MatrixXcd get_inv_Q(const Eigen::MatrixXcd &Q) override{
    return Q.adjoint();
  }
  void setup(Parameters &param){
  }
  EnvironmentProjector_Krylov(Parameters &param){
    setup(param);
  }
  virtual ~EnvironmentProjector_Krylov(){}
};

}//namespace

#endif
