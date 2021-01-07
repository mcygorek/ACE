#ifndef MODE_PROPAGATOR_DEFINED_H
#define MODE_PROPAGATOR_DEFINED_H

#include "FreePropagator.h"
#include "IF_OD_Dictionary.h"
#include "Smart_Ptr.h"

class ModePropagator: public FreePropagator{
public:
  int N_system;
//  int N_mode;
//  int factorization;

  std::vector<std::vector<Eigen::MatrixXcd> > A;
  std::vector<Eigen::MatrixXcd> env_ops;
  Eigen::MatrixXcd bath_init;

  int get_N_system()const{return N_system;}
  int get_N_mode()const{return bath_init.rows();}
  const Eigen::MatrixXcd &get_bath_init()const{return bath_init;}

  virtual void update(double t, double dt){
    FreePropagator::update(t,dt);

    int Ns=get_N_system();
    int Nm=get_N_mode();

    A.resize(Ns*Ns, std::vector<Eigen::MatrixXcd>(Ns*Ns, Eigen::MatrixXcd(Nm*Nm,Nm*Nm)));

    for(int nu1=0; nu1<Ns; nu1++){
      for(int mu1=0; mu1<Ns; mu1++){
        for(int nu2=0; nu2<Ns; nu2++){
          for(int mu2=0; mu2<Ns; mu2++){
            for(int xi1=0; xi1<Nm; xi1++){
              for(int xi_1=0; xi_1<Nm; xi_1++){
                for(int xi2=0; xi2<Nm; xi2++){
                  for(int xi_2=0; xi_2<Nm; xi_2++){
                    A[nu1*Ns+mu1][nu2*Ns+mu2](xi1*Nm+xi_1, xi2*Nm+xi_2)=
                    M( (nu1*Nm+xi1)*Nm*Ns+(mu1*Nm+xi_1),(nu2*Nm+xi2)*Nm*Ns+(mu2*Nm+xi_2));
                  }
                }
              }
            }
          }
        }
      }
    }

  }
  IF_OD_Dictionary get_dict(){
    return IF_OD_Dictionary(get_N_system());
  }

  ModePropagator(int Ns=2, 
                 const Eigen::MatrixXcd &binit=Eigen::MatrixXcd::Identity(2,2) )
    : N_system(Ns), bath_init(binit){
  }
  ModePropagator(int Ns, 
                 const Eigen::MatrixXcd &binit,
                 const Eigen::MatrixXcd &H, 
       std::vector<Eigen::MatrixXcd> env_mat=std::vector<Eigen::MatrixXcd>() ) 
    : N_system(Ns), bath_init(binit), env_ops(env_mat) {
    if(H.rows()!=Ns*get_N_mode() || H.cols()!=Ns*get_N_mode()){
      std::cerr<<"Error constructing ModePropagator: .rows()!=Ns*Nm || H.cols()!=Ns*Nm!"<<std::endl;
      exit(1);
    }
    set_Hamiltonian(H);
  }
  virtual ~ModePropagator(){}

};

typedef Smart_Ptr<ModePropagator> ModePropagatorPtr;

#endif
