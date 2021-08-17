#ifndef MODE_PROPAGATOR_DEFINED_H
#define MODE_PROPAGATOR_DEFINED_H

#include "FreePropagator.h"
#include "IF_OD_Dictionary.h"
#include "Smart_Ptr.h"
#include "Equilibrium.h"
#include "ReducedLiouvilleBasis.h"

class ModePropagator: public FreePropagator{
public:
  int N_system;
//  int N_mode;
//  int factorization;

  std::vector<std::vector<Eigen::MatrixXcd> > A;
  Eigen::MatrixXcd bath_init; 
  std::vector<Eigen::MatrixXcd> env_ops;
  Smart_Ptr<ReducedLiouvilleBasis> rBasis;


  struct low_pass_struct{
    bool use;
    double cutoff;
    double factor;
   
    double f(double w)const{
      double wabs=fabs(w);
      if(wabs<cutoff)return 0.;
      return factor*(wabs/cutoff-1.);
    }
    low_pass_struct():use(false), cutoff(1.){}
  }low_pass;

  struct continuum_subdiv_struct{
    int N;  //number of energy space subdivision
    Eigen::MatrixXcd dH; // dE of full interval times Operator Hdiag
    
    continuum_subdiv_struct() : N(0){}
  }continuum_subdiv;


  int get_N_system()const{return N_system;}
  int get_N_mode()const{return bath_init.rows();}
  const Eigen::MatrixXcd &get_bath_init()const{return bath_init;}


  void M_to_A(){
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

  virtual void update_single(double t, double dt)  { //override
    if(continuum_subdiv.N>0){
std::cout<<"TEST: subdiv: "<<continuum_subdiv.N<<std::endl;
std::cout<<continuum_subdiv.dH<<std::endl;
      int dimHSE=const_H.rows();
      if(continuum_subdiv.dH.rows()!=dimHSE){
        std::cerr<<"ModePropagator: continuum_subdiv.dH.rows()!=dimHSE!"<<std::endl;
        exit(1);
      }
      Eigen::MatrixXcd const_H_bck=const_H;
      M=Eigen::MatrixXcd::Zero(dimHSE*dimHSE, dimHSE*dimHSE);
      for(int n=0; n<continuum_subdiv.N; n++){
        double x=((double)n+0.5)/((double)continuum_subdiv.N);
        const_H=const_H_bck+continuum_subdiv.dH*(x-0.5);
        Eigen::MatrixXcd gen=Total_Generator(t,dt);
        M+=gen.exp();
      }
      const_H=const_H_bck;
      M=hs_rot.apply_Liouville(M/((double)continuum_subdiv.N));

    }else if(low_pass.use){
      std::cout<<"low pass cutoff: "<<low_pass.cutoff<<std::endl;
      Eigen::MatrixXcd gen=Total_Generator(t,dt);
      Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(gen);
//std::cout<<(solver.eigenvalues()/dt).transpose()<<std::endl;
      Eigen::VectorXcd diag=solver.eigenvalues();
      for(int i=0; i<diag.size(); i++){
        diag(i)=std::exp(diag(i)-low_pass.f(diag(i).imag()/dt)*dt);
      }
      M=solver.eigenvectors()*(diag.asDiagonal())*(solver.eigenvectors().adjoint());
      M=hs_rot.apply_Liouville(M);

    }else{
      Eigen::MatrixXcd gen=Total_Generator(t,dt);
      M=gen.exp();
      M=hs_rot.apply_Liouville(M);
    }
  }
  virtual void update(double t, double dt){
#ifdef DEBUG
std::cout<<"DEBUG: ModePropagator: update called: t="<<t<<std::endl;
#endif
    FreePropagator::update(t,dt);
    M_to_A();
  }

  IF_OD_Dictionary get_dict(){
    return IF_OD_Dictionary(get_N_system());
  }

  ModePropagator(int Ns=2, 
                 const Eigen::MatrixXcd &binit=Eigen::MatrixXcd::Identity(2,2) )
    : N_system(Ns), bath_init(binit){
    rBasis=new ReducedLiouvilleBasis();
  }
  ModePropagator(const FreePropagator &fprop,
      const Eigen::MatrixXcd &binit=Eigen::MatrixXcd::Identity(2,2),
      std::vector<Eigen::MatrixXcd> env_ops_=std::vector<Eigen::MatrixXcd>() )
    : FreePropagator(fprop), bath_init(binit), env_ops(env_ops_){

    N_system=fprop.get_dim()/bath_init.rows();
    if(N_system<2){
      std::cerr<<"ModePropagator(FreePropagator..): N_system<2!"<<std::endl;
      exit(1);
    }
    if(fprop.get_dim()!=N_system*bath_init.rows()){
      std::cerr<<"ModePropagator(FreePropagator..): fprop.get_dim()!=N_system*bath_init.rows()!"<<std::endl;
      exit(1);
    }
    rBasis=new ReducedLiouvilleBasis();
  }
  ModePropagator(int Ns, 
                 const Eigen::MatrixXcd &binit,
                 const Eigen::MatrixXcd &H, 
    std::vector<Eigen::MatrixXcd> env_mat=std::vector<Eigen::MatrixXcd>(), 
    const continuum_subdiv_struct cont_sub=continuum_subdiv_struct() ) 
    : N_system(Ns), bath_init(binit), env_ops(env_mat),
      continuum_subdiv(cont_sub) {

    if(H.rows()!=Ns*get_N_mode() || H.cols()!=Ns*get_N_mode()){
      std::cerr<<"Error constructing ModePropagator: .rows()!=Ns*Nm || H.cols()!=Ns*Nm!"<<std::endl;
      exit(1);
    }
    set_Hamiltonian(H);
    rBasis=new ReducedLiouvilleBasis();
  }
  virtual ~ModePropagator(){}

};

typedef Smart_Ptr<ModePropagator> ModePropagatorPtr;

#endif
