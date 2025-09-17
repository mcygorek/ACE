#include "ModePropagator.hpp"
#include "FreePropagator.hpp"
#include "IF_OD_Dictionary.hpp"
#include "Equilibrium.hpp"
#include "ReducedLiouvilleBasis.hpp"
#include "DummyException.hpp"


namespace ACE{

  void ModePropagator::M_to_A(){
    int Ns=get_N_system();
    int Nm=get_N_mode();

    A.resize(Ns*Ns, std::vector<Eigen::MatrixXcd>(Ns*Ns, Eigen::MatrixXcd(Nm*Nm,Nm*Nm)));

    if(fermion_sign_space<0){
      for(int nu1=0; nu1<Ns; nu1++){
        for(int mu1=0; mu1<Ns; mu1++){
          for(int nu2=0; nu2<Ns; nu2++){
            for(int mu2=0; mu2<Ns; mu2++){
              for(int xi2=0; xi2<Nm; xi2++){
                for(int xi_2=0; xi_2<Nm; xi_2++){
Eigen::Map<Eigen::MatrixXcd>(&A[nu1*Ns+mu1][nu2*Ns+mu2](0,xi2*Nm+xi_2),Nm,Nm).noalias()=  
Eigen::Map<Eigen::MatrixXcd, 0, Eigen::OuterStride<> >(
&M((nu1*Nm*Ns+mu1)*Nm, (nu2*Nm+xi2)*Nm*Ns+(mu2*Nm+xi_2)), Nm, Nm, 
Eigen::OuterStride<>(Nm*Ns));
/*
                  for(int xi1=0; xi1<Nm; xi1++){
                    for(int xi_1=0; xi_1<Nm; xi_1++){
A[nu1*Ns+mu1][nu2*Ns+mu2](xi1*Nm+xi_1, xi2*Nm+xi_2)=
   M( (nu1*Nm+xi1)*Nm*Ns+(mu1*Nm+xi_1),(nu2*Nm+xi2)*Nm*Ns+(mu2*Nm+xi_2));
                    }
                  }
*/
                }
              }
            }
          }
        }
      }
    }else{
      int mmod=Ns/pow(2,fermion_sign_space)/2;
      for(int nu1=0; nu1<Ns; nu1++){
        for(int mu1=0; mu1<Ns; mu1++){
          for(int nu2=0; nu2<Ns; nu2++){
            for(int mu2=0; mu2<Ns; mu2++){
              for(int xi2=0; xi2<Nm; xi2++){
                for(int xi_2=0; xi_2<Nm; xi_2++){
                  for(int xi1=0; xi1<Nm; xi1++){
                    for(int xi_1=0; xi_1<Nm; xi_1++){
                      double P=1.;
                      if((nu1/mmod)%2>0 && xi1%2>0)P*=-1.; 
                      if((mu1/mmod)%2>0 && xi_1%2>0)P*=-1.; 
A[nu1*Ns+mu1][nu2*Ns+mu2](xi1*Nm+xi_1, xi2*Nm+xi_2)=
   P*M( (nu1*Nm+xi1)*Nm*Ns+(mu1*Nm+xi_1),(nu2*Nm+xi2)*Nm*Ns+(mu2*Nm+xi_2));
                    }
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
    Eigen::MatrixXcd gen=Total_Generator(t,dt);
    M=gen.exp();
    M=hs_rot.apply_Liouville(M);
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
      EnvironmentOperators env_ops_)
    : FreePropagator(fprop), bath_init(binit), env_ops(env_ops_), 
      fermion_sign_space(-1){

    int N_system=fprop.get_dim()/bath_init.rows();
    if(N_system<2){
      std::cerr<<"ModePropagator(FreePropagator..): N_system<2!"<<std::endl;
      throw DummyException();
    }
    if(fprop.get_dim()!=N_system*bath_init.rows()){
      std::cerr<<"ModePropagator(FreePropagator..): fprop.get_dim()!=N_system*bath_init.rows()!"<<std::endl;
      throw DummyException();
    }
    rBasis=std::make_shared<ReducedLiouvilleBasis>();
  }

  ModePropagator::ModePropagator(
                 const Eigen::MatrixXcd &binit,
                 const Eigen::MatrixXcd &H, 
                 EnvironmentOperators env_mat)
    : bath_init(binit), env_ops(env_mat), fermion_sign_space(-1){

//    if(H.rows()!=Ns*get_N_mode() || H.cols()!=Ns*get_N_mode()){
//      std::cerr<<"Error constructing ModePropagator: .rows()!=Ns*Nm || H.cols()!=Ns*Nm!"<<std::endl;
//      throw DummyException();
//    }
    set_Hamiltonian(H);
    rBasis=std::make_shared<ReducedLiouvilleBasis>();
  }

}//namespace
