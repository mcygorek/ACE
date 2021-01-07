#ifndef FULL_PATH_SUM_MPS_DEFINED_H
#define FULL_PATH_SUM_MPS_DEFINED_H

#include "FullInfluenceFunctional_MPS.h"
#include "Propagator.h"
#include "Simulation_Results.h"

class FullPathSum_MPS{
public:

  Eigen::MatrixXcd rho;

  Eigen::MatrixXcd state; //laterally augmented density matrix


  Simulation_Results results;
  Output_Ops output_Op;

  void add_output_Op(const Eigen::MatrixXcd &op){
    output_Op.add(op);
  }


  void run(Propagator &prop, const FullInfluenceFunctional_MPS & FullIF,
    double ta, double dt, double te, const Eigen::MatrixXcd &initial_rho){

    if(initial_rho.rows()!=initial_rho.cols()){
      std::cerr<<"FullPathSum_MPS::calculate: initial_rho.rows()!=initial_rho.cols()!"<<std::endl;
      exit(1);
    }

    int N=initial_rho.rows();
    int NL=N*N;

    int n_max=(te-ta)/dt;
    if(FullIF.get_rank()<n_max+1){
      std::cerr<<"FullPathSum_MPS::calculate: Not enough steps calculated for influence functional!"<<std::endl;
      exit(1);
    }
 
    if(n_max<1){
      rho=initial_rho;
      return; 
    }

    if(FullIF.get_rank()<n_max+1){
      std::cerr<<"FullPathSum_MPS::calculate: FullIF.get_rank()<nmax+1!"<<std::endl;
      exit(1);
    }
    if(NL!=FullIF.get_dim(0)){
      std::cerr<<"FullPathSum_MPS::calculate: FullIF.get_dim(0)()!=NL!"<<std::endl;
      exit(1);
    }
    //initialize
    state.resize(NL, 1);
    for(int i=0; i<N; i++){
      for(int j=0; j<N; j++){
        state(i*N+j,0)=initial_rho(i,j);
      }
    }

    results.clear();
    results.resize(n_max+1);
    results.set(0, ta, output_Op, initial_rho);

    for(int n=0; n<n_max; n++){
      double t=ta+n*dt;
      prop.update(t,dt);

      const MPS_Matrix & a=FullIF.a[FullIF.a.size()-1-n];
      Eigen::MatrixXcd state2=Eigen::MatrixXcd::Zero(NL,a.dim_d1);
      std::cout<<"n: "<<n<<std::endl; //<<" from back: "<<FullIF.a.size()-1-n<<std::endl;
//      a.print_dims(); std::cout<<std::endl;
      for(int i=0; i<NL; i++){
        for(int j=0; j<NL; j++){
          for(int d1=0; d1<a.dim_d1; d1++){
            for(int d2=0; d2<a.dim_d2; d2++){
              state2(i,d1)+=prop.M(i,j)*a(j,d1,d2)*state(j,d2);
            }
          }
        }
      }
      state=state2;

      //extract reduced system density matrix
      rho=Eigen::MatrixXcd::Zero(N,N);
      const Eigen::MatrixXcd & c=FullIF.c[FullIF.a.size()-1-n-1];
      for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
          for(int d2=0; d2<c.cols(); d2++){
            rho(i,j) += c(i*N+j, d2) * state(i*N+j, d2);
          }
        }
      }
      results.set(n+1, ta+(n+1)*dt, output_Op, rho);
    }
  }


  void print_results(const std::string &fname)const{
    results.print(fname);
  }


  FullPathSum_MPS(Propagator &prop, const FullInfluenceFunctional_MPS & FullIF, 
    double ta, double dt, double te, const Eigen::MatrixXcd &initial_rho){
   
     run(prop, FullIF, ta, dt, te, initial_rho);
  }
  FullPathSum_MPS(){
  }
};

#endif
