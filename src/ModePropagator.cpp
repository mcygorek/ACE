#include "ModePropagator.hpp"
#include "FreePropagator.hpp"
#include "IF_OD_Dictionary.hpp"
#include "Equilibrium.hpp"
#include "ReducedLiouvilleBasis.hpp"


namespace ACE{

    double ModePropagator::low_pass_struct::f(double w)const{
      double wabs=fabs(w);
      if(wabs<cutoff)return 0.;
      return factor*(wabs/cutoff-1.);
    }


  void ModePropagator::M_to_A(){
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

  MPS_Matrix ModePropagator::get_MPS_Matrix()const{
    int N=get_N_system();
    int N_mode=get_N_mode();
    int NL=N*N;
    int ML=N_mode*N_mode;
    
    MPS_Matrix M(NL*NL, ML, ML);
    for(int i=0; i<NL; i++){
      for(int j=0; j<NL; j++){
        for(int d1=0; d1<ML; d1++){
          for(int d2=0; d2<ML; d2++){
            M(i*NL+j, d1, d2)=A[i][j](d2,d1);
          }
        }
      } 
    }
    return M;
  }

  void ModePropagator::update_single(double t, double dt)  { //override
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

  void ModePropagator::update(double t, double dt){
#ifdef DEBUG
std::cout<<"DEBUG: ModePropagator: update called: t="<<t<<std::endl;
#endif
    FreePropagator::update(t,dt);
    M_to_A();
  }

  ModePropagator::ModePropagator(const FreePropagator &fprop,
      const Eigen::MatrixXcd &binit,
      std::vector<Eigen::MatrixXcd> env_ops_)
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
    rBasis=std::make_shared<ReducedLiouvilleBasis>();
  }

  ModePropagator::ModePropagator(int Ns, 
                 const Eigen::MatrixXcd &binit,
                 const Eigen::MatrixXcd &H, 
    std::vector<Eigen::MatrixXcd> env_mat,
    const continuum_subdiv_struct cont_sub ) 
    : N_system(Ns), bath_init(binit), env_ops(env_mat),
      continuum_subdiv(cont_sub) {

    if(H.rows()!=Ns*get_N_mode() || H.cols()!=Ns*get_N_mode()){
      std::cerr<<"Error constructing ModePropagator: .rows()!=Ns*Nm || H.cols()!=Ns*Nm!"<<std::endl;
      exit(1);
    }
    set_Hamiltonian(H);
    rBasis=std::make_shared<ReducedLiouvilleBasis>();
  }

}//namespace
