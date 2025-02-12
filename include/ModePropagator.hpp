#ifndef MODE_PROPAGATOR_DEFINED_H
#define MODE_PROPAGATOR_DEFINED_H

#include "FreePropagator.hpp"
#include "IF_OD_Dictionary.hpp"
#include "Equilibrium.hpp"
#include "ReducedLiouvilleBasis.hpp"
#include "MPS_Matrix.hpp"
#include "EnvironmentOperators.hpp"
#include <memory>

namespace ACE{

class ModePropagator: public FreePropagator{
public:
  int N_system;
  int fermion_sign_space;

  std::vector<std::vector<Eigen::MatrixXcd> > A;
  Eigen::MatrixXcd bath_init; 
  EnvironmentOperators env_ops;
  std::shared_ptr<ReducedLiouvilleBasis> rBasis;



  inline int get_N_system()const{return N_system;}
  inline int get_N_mode()const{return bath_init.rows();}
  inline const Eigen::MatrixXcd &get_bath_init()const{return bath_init;}


  void M_to_A();

  virtual void update_single(double t, double dt); 
  virtual void update(double t, double dt);

  inline IF_OD_Dictionary get_dict(){
    return IF_OD_Dictionary(get_N_system());
  }
  
  MPS_Matrix get_MPS_Matrix()const;

  inline ModePropagator(int Ns=2, 
                 const Eigen::MatrixXcd &binit=Eigen::MatrixXcd::Identity(2,2) )
    : N_system(Ns), bath_init(binit), fermion_sign_space(-1){
    rBasis=std::make_shared<ReducedLiouvilleBasis>();
  }

  ModePropagator(const FreePropagator &fprop,
      const Eigen::MatrixXcd &binit=Eigen::MatrixXcd::Identity(2,2),
      EnvironmentOperators env_ops_=EnvironmentOperators() );
  
  ModePropagator(int Ns, 
                 const Eigen::MatrixXcd &binit,
                 const Eigen::MatrixXcd &H, 
                 EnvironmentOperators env_mat=EnvironmentOperators());
  
  virtual ~ModePropagator(){}
};

typedef std::shared_ptr<ModePropagator> ModePropagatorPtr;

}//namespace
#endif
