#ifndef FULL_PATH_SUM_DEFINED_H
#define FULL_PATH_SUM_DEFINED_H

#include "Tensor_Dense.h"
#include "Propagator.h"
#include "Simulation_Results.h"

class FullPathSum{
public:

  Eigen::MatrixXcd rho;

  Simulation_Results results;
  Output_Ops output_Op;

  void add_output_Op(const Eigen::MatrixXcd &op){
    output_Op.add(op);
  }


  void calculate( Propagator &prop, const Tensor & FullIF, 
      double ta, double dt, int n_max, const Eigen::MatrixXcd &initial_rho){


    if(initial_rho.rows()!=initial_rho.cols()){
      std::cerr<<"FullPathSum::calculate: initial_rho.rows()!=initial_rho.cols()!"<<std::endl;
      exit(1);
    }

    int N=initial_rho.rows();
    int NL=N*N;

    if(n_max<1){
      rho=initial_rho;
      return; 
    }

    if(FullIF.get_rank()<n_max+1){
      std::cerr<<"FullPathSum::calculate: FullIF.get_rank()<nmax+1!"<<std::endl;
      exit(1);
    }
    if(NL!=FullIF.get_dim(0)){
      std::cerr<<"FullPathSum::calculate: FullIF.get_dim(0)()!=NL!"<<std::endl;
      exit(1);
    }



    std::vector<Eigen::MatrixXcd> Mvec(n_max);
    for(size_t n=0; n<Mvec.size(); n++){
      double t=ta+n*dt;
      prop.update(t,dt);
      Mvec[n]=prop.M;
    }



    std::complex<double> *init_rho_vec=new std::complex<double>[NL];
    for(int i=0; i<N; i++){
      for(int j=0; j<N; j++){
        init_rho_vec[i*N+j]=initial_rho(i,j);
      }
    }

    std::complex<double> *rho_vec=new std::complex<double>[NL];
    for(int i=0; i<NL; i++)rho_vec[i]=0;

    int firsti=FullIF.get_rank()-(n_max+1);
    for(Tensor_Index ind(FullIF); !ind.done(FullIF); ind.increment_from_back(FullIF)){
      if(firsti>0 && ind[firsti-1]>0)break;
      std::complex<double> contrib=FullIF(ind);

      for(int r=firsti+1; r<ind.get_rank(); r++){
        contrib*=Mvec[ind.get_rank()-1-r](ind[r-1],ind[r]);
      }
      contrib*=init_rho_vec[ind.back()];
      rho_vec[ind[firsti]]+=contrib;
    }

    rho=Eigen::MatrixXcd::Zero(N,N);
    for(int i=0; i<N; i++){
      for(int j=0; j<N; j++){
        rho(i,j)=rho_vec[i*N+j];
      }
    }

    delete[] rho_vec;
    delete[] init_rho_vec;
}

  void run(Propagator &prop, const Tensor & FullIF,
    double ta, double dt, int n_max, const Eigen::MatrixXcd &initial_rho){
    results.clear();
    results.resize(n_max+1);

    results.set(0, ta, output_Op, initial_rho);
    for(size_t st=1; st<n_max+1; st++){
      calculate(prop, FullIF, ta, dt, st, initial_rho);
      results.set(st, ta+st*dt, output_Op, rho);
    }
  }

  void print_results(const std::string &fname)const{
    results.print(fname);
  }


  FullPathSum( Propagator &prop, const Tensor & FullIF, 
    double ta, double dt, int n_max, const Eigen::MatrixXcd &initial_rho){
   
     run(prop, FullIF, ta, dt, n_max, initial_rho);
  }
  FullPathSum(){
  }
};

#endif
