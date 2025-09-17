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
  int fermion_sign_space;

  std::vector<std::vector<Eigen::MatrixXcd> > A;
  Eigen::MatrixXcd bath_init; 
  EnvironmentOperators env_ops;
  std::shared_ptr<ReducedLiouvilleBasis> rBasis;



  inline int get_N_system()const{if(get_N_mode()<2)return 0; 
                                 return get_dim()/get_N_mode();}
  inline int get_N_mode()const{return bath_init.rows();}
  inline const Eigen::MatrixXcd &get_bath_init()const{return bath_init;}

  inline FreePropagator & get_fprop(){
    return (FreePropagator &)*this;
  }
  inline Eigen::MatrixXcd & get_initial(){
    return bath_init;
  }

  void M_to_A();

  virtual void update_single(double t, double dt); 
  virtual void update(double t, double dt);

  inline IF_OD_Dictionary get_dict(){
    return IF_OD_Dictionary(get_N_system());
  }
  
  MPS_Matrix get_MPS_Matrix()const;

  inline ModePropagator(int Ns=2, 
                 const Eigen::MatrixXcd &binit=Eigen::MatrixXcd::Identity(2,2) )
    : FreePropagator(Ns*binit.rows()), bath_init(binit), fermion_sign_space(-1){
    rBasis=std::make_shared<ReducedLiouvilleBasis>();
  }

  ModePropagator(const FreePropagator &fprop,
      const Eigen::MatrixXcd &binit=Eigen::MatrixXcd::Identity(2,2),
      EnvironmentOperators env_ops_=EnvironmentOperators() );
  
  ModePropagator(const Eigen::MatrixXcd &binit,
                 const Eigen::MatrixXcd &H, 
                 EnvironmentOperators env_mat=EnvironmentOperators());
  
  ModePropagator(const std::string & filename, const Eigen::MatrixXcd &initial)
    : FreePropagator(filename), bath_init(initial), fermion_sign_space(-1) {
    rBasis=std::make_shared<ReducedLiouvilleBasis>();
  }
  ModePropagator(Parameters & param, const Eigen::MatrixXcd &initial)
    : FreePropagator(param), bath_init(initial), fermion_sign_space(-1) {
    rBasis=std::make_shared<ReducedLiouvilleBasis>();
  }

  virtual ~ModePropagator(){}
};

typedef std::shared_ptr<ModePropagator> ModePropagatorPtr;

}//namespace
#endif
