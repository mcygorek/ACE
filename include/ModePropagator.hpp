#ifndef MODE_PROPAGATOR_DEFINED_H
#define MODE_PROPAGATOR_DEFINED_H

#include "FreePropagator.hpp"
#include "IF_OD_Dictionary.hpp"
#include "Equilibrium.hpp"
#include "ReducedLiouvilleBasis.hpp"
#include "MPS_Matrix.hpp"
#include <memory>

namespace ACE{

class ModePropagator: public FreePropagator{
public:
  int N_system;
//  int N_mode;
//  int factorization;

  std::vector<std::vector<Eigen::MatrixXcd> > A;
  Eigen::MatrixXcd bath_init; 
  std::vector<Eigen::MatrixXcd> env_ops;
  std::shared_ptr<ReducedLiouvilleBasis> rBasis;


  struct low_pass_struct{
    bool use;
    double cutoff;
    double factor;
   
    double f(double w)const;
    low_pass_struct():use(false), cutoff(1.){}
  }low_pass;


  struct continuum_subdiv_struct{
    int N;  //number of energy space subdivision
    Eigen::MatrixXcd dH; // dE of full interval times Operator Hdiag
    
    continuum_subdiv_struct() : N(0){}
  }continuum_subdiv;


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
    : N_system(Ns), bath_init(binit){
    rBasis=std::make_shared<ReducedLiouvilleBasis>();
  }

  ModePropagator(const FreePropagator &fprop,
      const Eigen::MatrixXcd &binit=Eigen::MatrixXcd::Identity(2,2),
      std::vector<Eigen::MatrixXcd> env_ops_=std::vector<Eigen::MatrixXcd>() );
  
  ModePropagator(int Ns, 
                 const Eigen::MatrixXcd &binit,
                 const Eigen::MatrixXcd &H, 
    std::vector<Eigen::MatrixXcd> env_mat=std::vector<Eigen::MatrixXcd>(), 
    const continuum_subdiv_struct cont_sub=continuum_subdiv_struct() );
  
  virtual ~ModePropagator(){}
};

typedef std::shared_ptr<ModePropagator> ModePropagatorPtr;

}//namespace
#endif
