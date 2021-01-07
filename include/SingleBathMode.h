#ifndef SINGLE_BATH_MODE_DEFINED_H
#define SINGLE_BATH_MODE_DEFINED_H

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>
#include "Constants.h"
#include "FreePropagator.h"
#include "OuterProduct.h"
#include <iostream>
#include <cstdlib>

/**
Structure to store the interaction between a N-level system with a 
single bath mode

*/

//Note: indices: H( nu1*dim_xi+ xi1, nu2*dim_xi+xi2)


class SingleBathMode{
public:
  bool ready;

  int N; //dimension of N-level system
  int M; //dimension of bath mode

  Eigen::MatrixXcd H; //Hamiltonian
  //Time evolution in Liouville space: (vector indices: alpha, tilde{alpha}
  std::vector<std::vector<Eigen::MatrixXcd> > A;
 
  int get_N()const{return N;}
  int get_M()const{return M;}

  void check_dimensions(){
    if(H.rows()!= get_N()*get_M() || H.cols()!=get_N()*get_M() ){
      std::cerr<<"SingleBathMode::check_dimensions: H: "<<H.rows()<<"x"<<H.cols()<<" N: "<<get_N()<<" M: "<<get_M()<<std::endl;
      exit(1);
    }
  }
  void calculateA_fullexp(double dt){
    check_dimensions();
    Eigen::MatrixXcd U=FreePropagator::Hamiltonian_to_Liouville_Generator(H, dt).exp();
    /*
      Note: One index of   has the format: 
          (nu*dim_xi+xi)*M*N+(mu*dim_xi+ xi_) 
             which we want to map to a  (alpha*dim_chi + chi)
    */
    int N=get_N(); int NL=N*N;
    int M=get_M(); int ML=M*M;
    A.resize(NL, std::vector<Eigen::MatrixXcd>(NL, Eigen::MatrixXcd(ML,ML)));

    for(int nu1=0; nu1<N; nu1++){
      for(int mu1=0; mu1<N; mu1++){
        for(int nu2=0; nu2<N; nu2++){
          for(int mu2=0; mu2<N; mu2++){
            for(int xi1=0; xi1<M; xi1++){
              for(int xi_1=0; xi_1<M; xi_1++){
                for(int xi2=0; xi2<M; xi2++){
                  for(int xi_2=0; xi_2<M; xi_2++){
                    A[nu1*N+mu1][nu2*N+mu2](xi1*M+xi_1, xi2*M+xi_2)= 
                    U( (nu1*M+xi1)*M*N+(mu1*M+xi_1),(nu2*M+xi2)*M*N+(mu2*M+xi_2));
                  }
                }
              }
            }
          }
        }
      }
    }
    ready=true;
  }

  void calculateA_factor1(double dt){
    check_dimensions();
    /*
      Note: One index of   has the format: 
          (nu*dim_xi+xi)*M*N+(mu*dim_xi+ xi_) 
             which we want to map to a  (alpha*dim_chi + chi)
    */
    std::complex<double> gamma(0.,-dt/Constants::hbar_in_meV_ps);
    int N=get_N(); int NL=N*N;
    int M=get_M(); int ML=M*M;
    A.resize(NL, std::vector<Eigen::MatrixXcd>(NL, Eigen::MatrixXcd(ML,ML)));

    //no-flip exponentials:
    std::vector<Eigen::MatrixXcd> Ud(N);
    for(int i=0; i<N; i++){
      Eigen::MatrixXcd Tmp(M,M);
      for(int xi1=0; xi1<M; xi1++){
        for(int xi2=0; xi2<M; xi2++){
          Tmp(xi1,xi2)=gamma*H(i*M+xi1, i*M+xi2);
        }
      }
      Ud[i]=Tmp.exp();
    }

    for(int nu1=0; nu1<N; nu1++){
      for(int mu1=0; mu1<N; mu1++){
        for(int nu2=0; nu2<N; nu2++){
          for(int mu2=0; mu2<N; mu2++){
            for(int xi1=0; xi1<M; xi1++){
              for(int xi_1=0; xi_1<M; xi_1++){
                for(int xi2=0; xi2<M; xi2++){
                  for(int xi_2=0; xi_2<M; xi_2++){
                    Eigen::MatrixXcd fac1;
                    if(nu1==nu2){
                      fac1=Ud[nu1];
                    }else{
                      fac1=Eigen::MatrixXcd::Zero(M,M);
                      for(int x=0; x<M; x++){
                        for(int y=0; y<M; y++){
                          for(int z=0; z<M; z++){
                            fac1(x,y)+=gamma*Ud[nu1](x,z)*H(nu1*M+z, nu2*M+y);
                          }
                        }
                      }
                    }
                    Eigen::MatrixXcd fac2;
                    if(mu1==mu2){
                      fac2=Ud[mu1].conjugate();
                    }else{
                      fac2=Eigen::MatrixXcd::Zero(M,M);
                      for(int x=0; x<M; x++){
                        for(int y=0; y<M; y++){
                          for(int z=0; z<M; z++){
                            fac2(x,y)+=-gamma*std::conj(Ud[mu1](x,z)*H(mu1*M+z, mu2*M+y));
                          }
                        }
                      }
                    }

                    A[nu1*N+mu1][nu2*N+mu2](xi1*M+xi_1, xi2*M+xi_2)= 
                      fac1(xi1,xi2)* fac2(xi_1, xi_2);
//                    U( (nu1*M+xi1)*M*N+(mu1*M+xi_1),(nu2*M+xi2)*M*N+(mu2*M+xi_2));
                  }
                }
              }
            }
          }
        }
      }
    }
    ready=true;
  }
  void calculateA(double dt, int factor=0){
    if(factor==1){
      calculateA_factor1(dt);  
    }else{
      calculateA_fullexp(dt);
    }
  }

  SingleBathMode(int N_=0, int M_=0, const Eigen::MatrixXcd &H_=Eigen::MatrixXcd())
   : N(N_), M(M_), H(H_){
    check_dimensions();
    ready=false;
  }

};



#endif
